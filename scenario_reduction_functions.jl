##### SCENARIO REDUCTION 
# This code uses timeseries of load as input 

using Pkg
# Pkg.add(["CSV", "DataFrames", "StatsBase", "Plots", "Distributions"])

using CSV
using DataFrames
using StatsBase
using Plots
using Plots.Measures
using LaTeXStrings


###############################################################################################################################

using CSV, DataFrames, StatsBase, Plots, Distributions

# -----------------------
# FUNCTIONs: Compute Scenarios, one for fluvius consumers and one for generated scenarios with scenario generator
# -----------------------
function compute_reduced_scenarios_fluvius_consumers(num_consumers::Int, num_buckets::Int)
    # --- Step 1: Read consumer IDs and select the first num_consumers ---
    consumer_ids_df = CSV.read("consumertype_IDs.csv", DataFrame)
    consumer_ids = consumer_ids_df[!, "type1_IDs"]
    selected_ids = consumer_ids[1:num_consumers]
    
    # --- Step 2: Define a function to load and aggregate a consumer load profile ---
    function load_consumer_profile(filename::String)
        df = CSV.read(filename, DataFrame, header=false)
        nrows, ncols = size(df)
        data_matrix = Array{Float64}(undef, nrows, ncols)
        for i in 1:nrows
            for j in 1:ncols
                data_matrix[i, j] = parse(Float64, string(df[i, j]))
            end
        end
        @assert size(data_matrix, 1) == 96 "Expected 96 rows per day in the consumer file."
        hourly_profile = Float64[]
        for day in 1:ncols
            day_values = data_matrix[:, day]
            day_matrix = reshape(day_values, 4, 24)   # 4 quarter-hours per hour, 24 hours
            hourly_values = sum(day_matrix, dims=1)
            append!(hourly_profile, vec(hourly_values))
        end
        return hourly_profile
    end

    # --- Step 3: Read the load profiles of the selected consumers ---
    consumer_loads = Matrix{Float64}(undef, 8760, num_consumers)
    load_dir = raw"C:\Users\milan\OneDrive - KU Leuven\2nd master\Thesis\2nd semester\Type1_consumer_Fluvius"
    for (i, id) in enumerate(selected_ids)
        filename = joinpath(load_dir, "Type1consumerLoad$(id).csv")
        println("Loading consumer ", id, " from file ", filename)
        consumer_loads[:, i] = load_consumer_profile(filename)
    end

    # --- Step 4: Build the timeseries matrix L ---
    T = 8760
    N = num_consumers
    L = Matrix{Float64}(undef, T, 3 + N)
    L[:, 1] = collect(1:T)                      # Column 1: hour index (1 to 8760)
    L[:, 2] = [mod(t - 1, 24) + 1 for t in 1:T]  # Column 2: hour of the day (1 to 24)
    L[:, 4:end] = consumer_loads                # Columns 4..(3+N): individual consumer load profiles
    L[:, 3] = sum(consumer_loads, dims=2)[:]      # Column 3: aggregated load (sum over all consumers)

    # Convert to DataFrame with clear column names
    col_names = vcat(["hour of the year", "hour of the day", "Aggregated load (kWh=kW)"], ["consumer_$(id)" for id in selected_ids])
    L_df = DataFrame(L, Symbol.(col_names))
    
    # Rename so that the column names match the scenario reduction code
    data = copy(L_df)
    data[!, "hour of the year"] = convert.(Int, data[!, "hour of the year"])
    data[!, "hour of the day"] = convert.(Int, data[!, "hour of the day"])
    
    # --- Scenario Reduction Code ---
    # Alias for column names
    hour = data[!, "hour of the year"]
    agg_load = data[!, "Aggregated load (kWh=kW)"]
    
    # Sort by aggregated load (for CDF)
    sorted_indices = sortperm(agg_load)
    sorted_agg_load = agg_load[sorted_indices]
    sorted_hour = hour[sorted_indices]
    
    # PDF via histogram
    hist = fit(Histogram, agg_load, nbins=100)
    pdf_values = hist.weights ./ sum(hist.weights)
    bin_centers = [mean(b) for b in hist.edges[1][1:end-1] .+ diff(hist.edges[1])/2]
    
    # CDF via empirical CDF
    ecdf_func = ecdf(agg_load)
    cdf_values = ecdf_func.(sorted_agg_load)
    
    # Define buckets with the input parameter num_buckets
    min_load = minimum(agg_load)
    max_load = maximum(agg_load)
    bucket_edges = range(min_load, stop=max_load, length=num_buckets+1)
    
    # Assign bucket to each scenario
    function assign_bucket(load)
        for i in 1:num_buckets
            if bucket_edges[i] ≤ load < bucket_edges[i+1] || (i == num_buckets && load == bucket_edges[i+1])
                return i
            end
        end
        return missing
    end
    data.bucket = assign_bucket.(agg_load)
    
    # Calculate bucket weights
    bucket_counts = combine(groupby(data, :bucket), nrow => :count)
    bucket_counts.bucket_weight = bucket_counts.count ./ sum(bucket_counts.count)
    
    # Cluster within buckets
    clustered_rows = DataFrame()
    for bucket_id in 1:num_buckets
        bucket_data = filter(row -> row.bucket == bucket_id, data)
        sorted_bucket = sort(bucket_data, Symbol("Aggregated load (kWh=kW)"))
        n = nrow(sorted_bucket)
        if n < 3
            # If there are too few elements, assign them all to the "mid" cluster.
            sorted_bucket.cluster_type = fill("mid", n)
            sorted_bucket.cluster_weight = ones(n)
        else
            # Determine the number of elements in the low cluster, at least 1 element
            n_low = max(round(Int, 0.10 * n), 1)           
            # Determine the number of elements for the mid cluster, at least 1 element
            n_mid = max(round(Int, 0.80 * n), 1)             
            # Ensure that at least 1 element remains for the high cluster
            if n_low + n_mid >= n
                n_mid = max(n - n_low - 1, 1)                # Adjusted: modification if low + mid >= n
            end
            low_idx = 1:n_low                               # Adjusted: low_idx based on n_low
            mid_idx = (n_low + 1):(n_low + n_mid)             # Adjusted: mid_idx based on n_low and n_mid
            high_idx = ((n_low + n_mid) + 1):n                # Adjusted: high_idx for the remaining rows

            sorted_bucket.cluster_type = vcat(
                fill("low", length(low_idx)),
                fill("mid", length(mid_idx)),
                fill("high", length(high_idx))
            )

            sorted_bucket.cluster_weight = similar(sorted_bucket.cluster_type, Float64)
            sorted_bucket.cluster_weight[low_idx] .= length(low_idx) / n
            sorted_bucket.cluster_weight[mid_idx] .= length(mid_idx) / n
            sorted_bucket.cluster_weight[high_idx] .= length(high_idx) / n
        end
        bw = bucket_counts.bucket_weight[bucket_counts.bucket .== bucket_id]
        if isempty(bw)
            bucket_w = 0.0
        else
            bucket_w = bw[1]
        end

        sorted_bucket.bucket_weight = fill(bucket_w, n)
        sorted_bucket.final_weight = sorted_bucket.bucket_weight .* sorted_bucket.cluster_weight
        sorted_bucket.centroid = falses(n)
        grouped = groupby(sorted_bucket, :cluster_type)
        for g in grouped
            μ = mean(g[!, Symbol("Aggregated load (kWh=kW)")])
            idx = findmin(abs.(g[!, Symbol("Aggregated load (kWh=kW)")] .- μ))[2]
            hour_val = g[idx, Symbol("hour of the year")]
            clustertype_val = g[idx, :cluster_type]
            row_idx = findfirst(r -> r[Symbol("hour of the year")] == hour_val && r.cluster_type == clustertype_val, eachrow(sorted_bucket))
            if !isnothing(row_idx)
                sorted_bucket[row_idx, :centroid] = true
            end
        end
        append!(clustered_rows, sorted_bucket)
    end
    
    # Select representative scenarios (centroids)
    representative_scenarios = filter(:centroid => ==(true), clustered_rows)
    output_df = select(representative_scenarios,
        :"hour of the year", "Aggregated load (kWh=kW)",
        :bucket, :cluster_type, :bucket_weight, :cluster_weight, :final_weight)
    
    # Return the original aggregated load vector (for metrics) and the output_df
    return (agg_load=agg_load, output_df=output_df)
end
function compute_reduced_scenarios_from_matrix(
    consumer_loads::Matrix{Float64},
    selected_ids::Vector{Int},
    num_buckets::Int)
    num_consumers = size(consumer_loads, 2)
    T = 8760

    # Step 1: Build the time series matrix L
    L = Matrix{Float64}(undef, T, 3 + num_consumers)
    L[:, 1] = collect(1:T)                           # Column 1: hour of the year
    L[:, 2] = [mod(t - 1, 24) + 1 for t in 1:T]      # Column 2: hour of the day
    L[:, 4:end] = consumer_loads                    # Columns 4..: individual consumer profiles
    L[:, 3] = sum(consumer_loads, dims=2)[:]        # Column 3: aggregated load

    # Step 2: Convert to DataFrame
    col_names = vcat(["hour of the year", "hour of the day", "Aggregated load (kWh=kW)"], 
                     ["consumer_$(id)" for id in selected_ids])
    L_df = DataFrame(L, Symbol.(col_names))

    # Step 3: Rename and prepare
    data = copy(L_df)
    data[!, "hour of the year"] = convert.(Int, data[!, "hour of the year"])
    data[!, "hour of the day"] = convert.(Int, data[!, "hour of the day"])
    agg_load = data[!, "Aggregated load (kWh=kW)"]

    # Step 4: Bucket assignment
    min_load = minimum(agg_load)
    max_load = maximum(agg_load)
    bucket_edges = range(min_load, stop=max_load, length=num_buckets + 1)

    function assign_bucket(load)
        for i in 1:num_buckets
            if bucket_edges[i] ≤ load < bucket_edges[i+1] || (i == num_buckets && load == bucket_edges[i+1])
                return i
            end
        end
        return missing
    end
    data.bucket = assign_bucket.(agg_load)

    # Step 5: Compute weights
    bucket_counts = combine(groupby(data, :bucket), nrow => :count)
    bucket_counts.bucket_weight = bucket_counts.count ./ sum(bucket_counts.count)

    # Step 6: Clustering within each bucket
    clustered_rows = DataFrame()
    for bucket_id in 1:num_buckets
        bucket_data = filter(row -> row.bucket == bucket_id, data)
        sorted_bucket = sort(bucket_data, Symbol("Aggregated load (kWh=kW)"))
        n = nrow(sorted_bucket)
        if n < 3
            sorted_bucket.cluster_type = fill("mid", n)
            sorted_bucket.cluster_weight = ones(n)
        else
            n_low = max(round(Int, 0.10 * n), 1)
            n_mid = max(round(Int, 0.80 * n), 1)
            if n_low + n_mid >= n
                n_mid = max(n - n_low - 1, 1)
            end
            low_idx = 1:n_low
            mid_idx = (n_low + 1):(n_low + n_mid)
            high_idx = ((n_low + n_mid) + 1):n

            sorted_bucket.cluster_type = vcat(
                fill("low", length(low_idx)),
                fill("mid", length(mid_idx)),
                fill("high", length(high_idx))
            )

            sorted_bucket.cluster_weight = similar(sorted_bucket.cluster_type, Float64)
            sorted_bucket.cluster_weight[low_idx] .= length(low_idx) / n
            sorted_bucket.cluster_weight[mid_idx] .= length(mid_idx) / n
            sorted_bucket.cluster_weight[high_idx] .= length(high_idx) / n
        end

        bw = bucket_counts.bucket_weight[bucket_counts.bucket .== bucket_id]
        bucket_w = isempty(bw) ? 0.0 : bw[1]

        sorted_bucket.bucket_weight = fill(bucket_w, n)
        sorted_bucket.final_weight = sorted_bucket.bucket_weight .* sorted_bucket.cluster_weight
        sorted_bucket.centroid = falses(n)

        grouped = groupby(sorted_bucket, :cluster_type)
        for g in grouped
            μ = mean(g[!, Symbol("Aggregated load (kWh=kW)")])
            idx = findmin(abs.(g[!, Symbol("Aggregated load (kWh=kW)")] .- μ))[2]
            hour_val = g[idx, Symbol("hour of the year")]
            clustertype_val = g[idx, :cluster_type]
            row_idx = findfirst(r -> r[Symbol("hour of the year")] == hour_val && 
                                     r.cluster_type == clustertype_val, eachrow(sorted_bucket))
            if !isnothing(row_idx)
                sorted_bucket[row_idx, :centroid] = true
            end
        end

        append!(clustered_rows, sorted_bucket)
    end

    # Step 7: Select centroids
    representative_scenarios = filter(:centroid => ==(true), clustered_rows)
    output_df = select(representative_scenarios,
        :"hour of the year", "Aggregated load (kWh=kW)",
        :bucket, :cluster_type, :bucket_weight, :cluster_weight, :final_weight)

    return (agg_load=agg_load, output_df=output_df)
end
function compute_reduced_scenarios_of_generated_scenarios(num_buckets::Int,generated_active_load_matrix::Matrix{Float64},generated_reactive_load_matrix::Matrix{Float64})

    # Build the timeseries matrix L ---
    T = size(generated_active_load_matrix, 1)
    N = size(generated_active_load_matrix, 2)
    L = Matrix{Float64}(undef, T, 3 + N)
    bus_id_range=1:N
    L[:, 1] = collect(1:T)                      # Column 1: hour index (1 to 8760)
    L[:, 2] = [mod(t - 1, 24) + 1 for t in 1:T]  # Column 2: hour of the day (1 to 24)
    L[:, 4:end] = generated_active_load_matrix               # Columns 4..(3+N): individual consumer load profiles
    L[:, 3] = sum(generated_active_load_matrix, dims=2)[:]      # Column 3: aggregated load (sum over all consumers)

    # Convert to DataFrame with clear column names
    col_names = vcat(["hour of the year", "hour of the day", "Aggregated load (MWh=MW)"], ["bus/consumer_$(id)" for id in bus_id_range])
    L_df = DataFrame(L, Symbol.(col_names))
    
    # transform first two columns into integers
    data = copy(L_df)
    data[!, "hour of the year"] = convert.(Int, data[!, "hour of the year"])
    data[!, "hour of the day"] = convert.(Int, data[!, "hour of the day"])
    
    # --- Scenario Reduction Code ---
    # Alias for column names
    hour = data[!, "hour of the year"]
    agg_load = data[!, "Aggregated load (MWh=MW)"]
    
    # Sort by aggregated load (for CDF)
    sorted_indices = sortperm(agg_load)
    sorted_agg_load = agg_load[sorted_indices]
    sorted_hour = hour[sorted_indices]
    
    # PDF via histogram
    hist = fit(Histogram, agg_load, nbins=100)
    pdf_values = hist.weights ./ sum(hist.weights)
    bin_centers = [mean(b) for b in hist.edges[1][1:end-1] .+ diff(hist.edges[1])/2]
    
    # CDF via empirical CDF
    ecdf_func = ecdf(agg_load)
    cdf_values = ecdf_func.(sorted_agg_load)
    
    # Define buckets with the input parameter num_buckets
    min_load = minimum(agg_load)
    max_load = maximum(agg_load)
    bucket_edges = range(min_load, stop=max_load, length=num_buckets+1)
    
    # Assign bucket to each scenario
    function assign_bucket(load)
        for i in 1:num_buckets
            if bucket_edges[i] ≤ load < bucket_edges[i+1] || (i == num_buckets && load == bucket_edges[i+1])
                return i
            end
        end
        return missing
    end
    data.bucket = assign_bucket.(agg_load)
    
    # Calculate bucket weights
    bucket_counts = combine(groupby(data, :bucket), nrow => :count)
    bucket_counts.bucket_weight = bucket_counts.count ./ sum(bucket_counts.count)
    
    # Cluster within buckets
    clustered_rows = DataFrame()
    for bucket_id in 1:num_buckets
        bucket_data = filter(row -> row.bucket == bucket_id, data)
        sorted_bucket = sort(bucket_data, Symbol("Aggregated load (MWh=MW)"))
        n = nrow(sorted_bucket)
        if n < 3
            # If there are too few elements, assign them all to the "mid" cluster.
            sorted_bucket.cluster_type = fill("mid", n)
            sorted_bucket.cluster_weight = ones(n)
        else
            # Determine the number of elements in the low cluster, at least 1 element
            n_low = max(round(Int, 0.10 * n), 1)           # Adjusted: use of max(..., 1)
            # Determine the number of elements for the mid cluster, at least 1 element
            n_mid = max(round(Int, 0.80 * n), 1)             # Adjusted: use of max(..., 1)
            # Ensure that at least 1 element remains for the high cluster
            if n_low + n_mid >= n
                n_mid = max(n - n_low - 1, 1)                # Adjusted: modification if low + mid >= n
            end
            low_idx = 1:n_low                               # Adjusted: low_idx based on n_low
            mid_idx = (n_low + 1):(n_low + n_mid)             # Adjusted: mid_idx based on n_low and n_mid
            high_idx = ((n_low + n_mid) + 1):n                # Adjusted: high_idx for the remaining rows

            sorted_bucket.cluster_type = vcat(
                fill("low", length(low_idx)),
                fill("mid", length(mid_idx)),
                fill("high", length(high_idx))
            )

            sorted_bucket.cluster_weight = similar(sorted_bucket.cluster_type, Float64)
            sorted_bucket.cluster_weight[low_idx] .= length(low_idx) / n
            sorted_bucket.cluster_weight[mid_idx] .= length(mid_idx) / n
            sorted_bucket.cluster_weight[high_idx] .= length(high_idx) / n
        end
        bw = bucket_counts.bucket_weight[bucket_counts.bucket .== bucket_id]
        if isempty(bw)
            bucket_w = 0.0
        else
            bucket_w = bw[1]
        end

        sorted_bucket.bucket_weight = fill(bucket_w, n)
        sorted_bucket.final_weight = sorted_bucket.bucket_weight .* sorted_bucket.cluster_weight
        sorted_bucket.centroid = falses(n)
        grouped = groupby(sorted_bucket, :cluster_type)
        for g in grouped
            μ = mean(g[!, Symbol("Aggregated load (MWh=MW)")])
            idx = findmin(abs.(g[!, Symbol("Aggregated load (MWh=MW)")] .- μ))[2]
            hour_val = g[idx, Symbol("hour of the year")]
            clustertype_val = g[idx, :cluster_type]
            row_idx = findfirst(r -> r[Symbol("hour of the year")] == hour_val && r.cluster_type == clustertype_val, eachrow(sorted_bucket))
            if !isnothing(row_idx)
                sorted_bucket[row_idx, :centroid] = true
            end
        end
        append!(clustered_rows, sorted_bucket)
    end
    
    # Select representative scenarios (centroids)
    representative_scenarios = filter(:centroid => ==(true), clustered_rows)
    output_df = select(representative_scenarios,
        :"hour of the year", "Aggregated load (MWh=MW)",
        :bucket, :cluster_type, :bucket_weight, :cluster_weight, :final_weight)
    
    # Return the original aggregated load vector (for metrics) and the output_df
    return (agg_load=agg_load, output_df=output_df)
end
function load_consumer_profile(filename::String)
    df = CSV.read(filename, DataFrame, header=false)
    nrows, ncols = size(df)
    data_matrix = Array{Float64}(undef, nrows, ncols)
    for i in 1:nrows, j in 1:ncols
        data_matrix[i, j] = parse(Float64, string(df[i, j]))
    end
    @assert size(data_matrix, 1) == 96 "Expected 96 rows per day in the consumer file."
    hourly_profile = Float64[]
    for day in 1:ncols
        day_values = data_matrix[:, day]
        day_matrix = reshape(day_values, 4, 24)
        hourly_values = sum(day_matrix, dims=1)
        append!(hourly_profile, vec(hourly_values))
    end
    return hourly_profile
end
function casestudy4_reducedscenarios()
    csv_path = "fluvius load profile (consumer1) and POLA network_timeseries_matrix.csv"
    data = CSV.read(csv_path, DataFrame)

    # Replace comma with dot and parse as Float64
    data[!, "Aggregated load (kWh=kW)"] = parse.(Float64, replace.(data[!, "Aggregated load (kWh=kW)"], "," => "."))

    # Apply same conversion to other relevant columns
    for col in names(data)[4:end]
        data[!, col] = parse.(Float64, replace.(data[!, col], "," => "."))
    end

    # Alias for column names
    hour = data[!, "hour of the year"]
    agg_load = data[!, "Aggregated load (kWh=kW)"]

    # Sort by aggregated load for CDF
    sorted_indices = sortperm(agg_load)
    sorted_agg_load = agg_load[sorted_indices]
    sorted_hour = hour[sorted_indices]

    # CDF via empirical CDF
    ecdf_func = ecdf(agg_load)
    cdf_values = ecdf_func.(sorted_agg_load)


    # Step 1: Define buckets
    num_buckets = 5
    min_load = minimum(agg_load)
    max_load = maximum(agg_load)
    bucket_edges = range(min_load, stop=max_load, length=num_buckets+1)


    # Assign bucket per scenario
    function assign_bucket(load)
        for i in 1:num_buckets
            if bucket_edges[i] ≤ load < bucket_edges[i+1] || (i == num_buckets && load == bucket_edges[i+1])
                return i
            end
        end
        return missing
    end

    data.bucket = assign_bucket.(agg_load)

    # Compute bucket weights
    bucket_counts = combine(groupby(data, :bucket), nrow => :count)
    bucket_counts.bucket_weight = bucket_counts.count ./ sum(bucket_counts.count)

    # Step 2: Cluster each bucket into 3 parts (low, mid, high) and determine centroids
    clustered_rows = DataFrame()

    for bucket_id in 1:length(bucket_counts.bucket)
        bucket_data = filter(row -> row.bucket == bucket_id, data)
        sorted_bucket = sort(bucket_data, Symbol("Aggregated load (kWh=kW)"))
        n = nrow(sorted_bucket)
        if n < 3
            # If there are too few elements, assign them all to the "mid" cluster.
            sorted_bucket.cluster_type = fill("mid", n)
            sorted_bucket.cluster_weight = ones(n)
        else
            # Determine the number of elements in the low cluster, at least 1 element
            n_low = max(round(Int, 0.10 * n), 1)           # Adjusted: use of max(..., 1)
            # Determine the number of elements for the mid cluster, at least 1 element
            n_mid = max(round(Int, 0.80 * n), 1)             # Adjusted: use of max(..., 1)
            # Ensure that at least 1 element remains for the high cluster
            if n_low + n_mid >= n
                n_mid = max(n - n_low - 1, 1)                # Adjusted: modification if low + mid >= n
            end
            low_idx = 1:n_low                               # Adjusted: low_idx based on n_low
            mid_idx = (n_low + 1):(n_low + n_mid)             # Adjusted: mid_idx based on n_low and n_mid
            high_idx = ((n_low + n_mid) + 1):n                # Adjusted: high_idx for the remaining rows

            sorted_bucket.cluster_type = vcat(
                fill("low", length(low_idx)),
                fill("mid", length(mid_idx)),
                fill("high", length(high_idx))
            )

            sorted_bucket.cluster_weight = similar(sorted_bucket.cluster_type, Float64)
            sorted_bucket.cluster_weight[low_idx] .= length(low_idx) / n
            sorted_bucket.cluster_weight[mid_idx] .= length(mid_idx) / n
            sorted_bucket.cluster_weight[high_idx] .= length(high_idx) / n
        end
        bw = bucket_counts.bucket_weight[bucket_counts.bucket .== bucket_id]
        if isempty(bw)
            bucket_w = 0.0
        else
            bucket_w = bw[1]
        end

        sorted_bucket.bucket_weight = fill(bucket_w, n)
        sorted_bucket.final_weight = sorted_bucket.bucket_weight .* sorted_bucket.cluster_weight
        sorted_bucket.centroid = falses(n)
        grouped = groupby(sorted_bucket, :cluster_type)
        for g in grouped
            μ = mean(g[!, Symbol("Aggregated load (kWh=kW)")])
            idx = findmin(abs.(g[!, Symbol("Aggregated load (kWh=kW)")] .- μ))[2]
            hour_val = g[idx, Symbol("hour of the year")]
            clustertype_val = g[idx, :cluster_type]
            row_idx = findfirst(r -> r[Symbol("hour of the year")] == hour_val && r.cluster_type == clustertype_val, eachrow(sorted_bucket))
            if !isnothing(row_idx)
                sorted_bucket[row_idx, :centroid] = true
            end
        end
        append!(clustered_rows, sorted_bucket)
    end


    # Select 15 representative scenarios (centroids)
    representative_scenarios = filter(:centroid => ==(true), clustered_rows)

    # Keep relevant columns
    output_df = select(representative_scenarios,
        :"hour of the year",
        "Aggregated load (kWh=kW)",
        :bucket,
        :cluster_type,
        :bucket_weight,
        :cluster_weight,
        :final_weight
    )
    println("output_df: ", output_df)

    # Show sum of final weights
    println("Sum of the 15 scenario weights: ", sum(output_df.final_weight))

    # Sort and export
    sort!(output_df, [:bucket, :cluster_type])
    return output_df
end
# -----------------------
# FUNCTION: Compute Metrics
# -----------------------
function compute_metrics(agg_load::Vector{Float64}, output_df::DataFrame)
    T = length(agg_load)

    # Annual Consumption
    annual_orig = sum(agg_load)
    annual_scen = 8760 * sum(output_df[!, "Aggregated load (kWh=kW)"] .* output_df.final_weight)
    annual_error = abs(annual_orig - annual_scen) / annual_orig

    # Peak Consumption
    peak_orig = maximum(agg_load)
    min_orig = minimum(agg_load)
    peak_scen = maximum(output_df[!, "Aggregated load (kWh=kW)"])
    peak_error = abs(peak_orig - peak_scen) / abs(peak_orig-min_orig)

    # Minimum (lowest) Load
    min_orig = minimum(agg_load)
    min_scen = minimum(output_df[!, "Aggregated load (kWh=kW)"])
    min_error = abs(min_orig - min_scen) / abs(peak_orig-min_orig)

    # Average Load
    avg_orig = mean(agg_load)
    avg_scen = sum(output_df[!, "Aggregated load (kWh=kW)"] .* output_df.final_weight)
    avg_error = abs(avg_orig - avg_scen) / avg_orig

    return (
        annual_error=annual_error,
        peak_error=peak_error,
        min_error=min_error,
        avg_error=avg_error,
        annual_orig=annual_orig,
        annual_scen=annual_scen,
        peak_orig=peak_orig,
        peak_scen=peak_scen,
        min_orig=min_orig,
        min_scen=min_scen,
        avg_orig=avg_orig,
        avg_scen=avg_scen
    )
end




#FUNCTION: fanchart
"""
    fanchart!(x, quantiles_matrix; label="", color=RGBA(0,0,1,0.3))

Draws a fanchart with 5 quantile levels: q10, q25, q50 (median), q75, q90.
"""

function fanchart!(plt, x::AbstractVector, y::AbstractMatrix; label="", color=RGBA(0, 0, 1, 0.5))
    @assert size(y, 2) == 5 "Matrix y must have 5 columns for quantiles: q10, q25, q50, q75, q90"
    
    plot!(plt, x, y[:, 3], label=label, color=color, linewidth=2)  # Median (q50)
    plot!(plt, x, y[:, 2:4], fillrange=y[:, 4], fillalpha=0.3, color=color, label="")  # IQR band
    plot!(plt, x, y[:, 1:5:5], fillrange=y[:, 5], fillalpha=0.15, color=color, label="")  # 10–90% band
end