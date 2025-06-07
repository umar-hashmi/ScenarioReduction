#this script consists of the code used for the case studies in the paper
using CSV, DataFrames, StatsBase, Plots, Distributions, Random
using PowerModels, JuMP, Ipopt

# FUNCTIONs (these two files need to be in the same folder as this script): 
include("scenario_reduction_functions.jl") #contribution of this paper
include("load_scenario_generation.jl") #load scenario generation functions based on a multivariate Gaussian distribution with Cholesky decomposition, as introduced in https://doi.org/10.1016/j.segan.2023.101069
include("CreatePMDDictionary.jl") #create PMD dictionary functions, author: Alexander Hoogsteyn
###########################################################################################################
###########################################################################################################
#Case study 1: Number of buckets
###########################################################################################################
###########################################################################################################

#Load duration curve (figure 2 in the paper)

###########################################################################################################
num_consumers_ldc = 100 # 100 consumers get selected randomly (first 100 from the dataset)
num_buckets_ldc = 5
res_ldc = compute_reduced_scenarios_fluvius_consumers(num_consumers_ldc, num_buckets_ldc)
agg_load_ldc = res_ldc.agg_load
output_df_ldc = res_ldc.output_df
output_df_ldc[!, "hour of the day"] = mod.(output_df_ldc[!, "hour of the year"] .- 1, 24) .+ 1
cols = names(output_df_ldc)
reordered = vcat(cols[1], "hour of the day", setdiff(cols, [cols[1], "hour of the day"]))
output_df_ldc = output_df_ldc[:, reordered]

println("\n--- Results for num_consumers = $num_consumers_ldc, num_buckets = $num_buckets_ldc ---")
println(output_df_ldc)
# Original Load Duration Curve:
# Sort the aggregated load in descending order
ldc_orig = sort(agg_load_ldc, rev=true)


# Reduced Load Duration Curve:
# Use output_df_ldc, which contains the representative scenarios
# Calculate for each representative scenario the number of hours it represents:
hour_counts = 8760 .* output_df_ldc.final_weight

# Sort the representative scenarios by their aggregated load (in descending order)
reduced_sorted = sort(output_df_ldc, "Aggregated load (kWh=kW)", rev=true)
# Calculate the cumulative hours based on the final_weight, so that a step-plot can be made:
cum_hours = cumsum(8760 .* reduced_sorted.final_weight)

# Create a step curve based on the representative load values.
# Since the reduced set contains discrete points, we plot the load value over the cumulative number of hours
# === Start building the plot ===
plt = Plots.plot(ldc_orig,
     xlabel="Hours",
     ylabel="Aggregated Load (kW)",
     #title="Load Duration Curve",
     label="Original LDC",
     linewidth=2)

# Add the reduced (step-wise) LDC
Plots.plot!(plt, cum_hours, reduced_sorted[!, "Aggregated load (kWh=kW)"],
      seriestype=:step,
      color=:green,
      linewidth=3,
      label="Reduced LDC")

# Add star markers on the reduced LDC
scatter!(plt, cum_hours, reduced_sorted[!, "Aggregated load (kWh=kW)"],
         marker=:star5,
         color=:yellow,
         markerstrokecolor=:black,
         markerstrokewidth=2,
         markersize=10,
         label="Representative scenarios")

# Annotate each star with its scenario number (1 to 15)
# Annotate each star with a number, connected by a line
for (i, (x, y)) in enumerate(zip(cum_hours, reduced_sorted[!, "Aggregated load (kWh=kW)"]))
    # bigger offset to avoid overlap
    x_offset = x + 250
    y_offset = y + 10

    # Dotted line from star to number
    plot!([x, x_offset], [y, y_offset], color=:gray, linestyle=:dot, label="")

    # Scenario number in red
    annotate!(x_offset, y_offset, Plots.text(string(length(cum_hours) - i + 1), :red, 12))
end


# Show the final plot
display(plt)

###########################################################################################################

#Fanchart of peak errors for MC of N_L=50 (figure 3 in the paper)

###########################################################################################################
# Read and store all 300 profiles in memory
consumer_ids_df = CSV.read("consumertype_IDs.csv", DataFrame; delim=';')
all_ids = Int64.(consumer_ids_df[1:300, "type1_IDs"])
load_dir = raw"C:\Users\milan\OneDrive - KU Leuven\2nd master\Thesis\2nd semester\Type1_consumer_Fluvius"
all_profiles = Dict{Int, Vector{Float64}}()

for id in all_ids
    filename = joinpath(load_dir, "Type1consumerLoad$(id).csv")
    all_profiles[id] = load_consumer_profile(filename)
end

# === Monte Carlo settings ===
num_trials = 100
bucket_range = 1:15
consumer_sizes = [15, 50, 100, 200]
results_mc = Dict{Int, DataFrame}()
# === Run Monte Carlo ===
# === Run Monte Carlo ===
for N in consumer_sizes
    println("Running Monte Carlo for N = $N")
    df = DataFrame(bucket=Int[], peak_error=Float64[], min_error=Float64[], avg_error=Float64[], annual_error=Float64[], trial=Int[])
    for trial in 1:num_trials
        println("Trial $trial")
        selected_ids = sample(all_ids, N, replace=false)
        selected_profiles = [all_profiles[id] for id in selected_ids]
        consumer_loads = hcat(selected_profiles...)  # 8760 x N

        for B in bucket_range
            res = compute_reduced_scenarios_from_matrix(consumer_loads, selected_ids, B)
            metrics = compute_metrics(res.agg_load, res.output_df)
            push!(df, (B, metrics.peak_error, metrics.min_error, metrics.avg_error, metrics.annual_error, trial)) 
        end
    end
    results_mc[N] = df
end


###plotting fancharts for different number of consumers
# === Plotting ===

peak_colors = Dict(
    15 => RGBA(0.9, 0.1, 0.1, 0.5),   # Red
    50 => RGBA(1.0, 0.5, 0.0, 0.5),   # Orange
    100 => RGBA(0.8, 0.6, 0.0, 0.5),  # yellow
    200 => RGBA(0.6, 0.8, 0.0, 0.5)   # green
)

min_colors = Dict(
    15 => RGBA(0.2, 0.2, 1.0, 0.5),   # blue
    50 => RGBA(0.4, 0.6, 1.0, 0.5),   # light blue
    100 => RGBA(0.6, 0.8, 1.0, 0.5),  # Sky blue
    200 => RGBA(0.7, 0.9, 1.0, 0.5)   # lightest blue
)

avg_colors = Dict(
    15 => RGBA(0.5, 0.0, 0.5, 0.5),   # Purple
    50 => RGBA(0.7, 0.3, 0.7, 0.5),   # Light Purple
    100 => RGBA(0.8, 0.4, 0.8, 0.5),  # Pinkish
    200 => RGBA(0.9, 0.5, 0.9, 0.5)   # Light Pink
)

plt = plot(xlabel="Number of Buckets B", ylabel="Error")
N = 50
# Plotting for N = 50
df = results_mc[N]

# Peak error
grouped_peak = combine(groupby(df, :bucket)) do sub
    q = quantile(sub.peak_error, [0.1, 0.25, 0.5, 0.75, 0.9])
    (; bucket=sub.bucket[1], q10=q[1], q25=q[2], q50=q[3], q75=q[4], q90=q[5])
end
y_peak = Matrix(grouped_peak[:, [:q10, :q25, :q50, :q75, :q90]])
fanchart!(plt, grouped_peak.bucket, y_peak,
          label="Upper Peak Error",color=:red)

# Min error
grouped_min = combine(groupby(df, :bucket)) do sub
    q = quantile(sub.min_error, [0.1, 0.25, 0.5, 0.75, 0.9])
    (; bucket=sub.bucket[1], q10=q[1], q25=q[2], q50=q[3], q75=q[4], q90=q[5])
end
y_min = Matrix(grouped_min[:, [:q10, :q25, :q50, :q75, :q90]])
fanchart!(plt, grouped_min.bucket, y_min,
          label="Lower Peak Error",color=:blue)

# Average error
grouped_avg = combine(groupby(df, :bucket)) do sub
    q = quantile(sub.avg_error, [0.1, 0.25, 0.5, 0.75, 0.9])
    (; bucket=sub.bucket[1], q10=q[1], q25=q[2], q50=q[3], q75=q[4], q90=q[5])
end
y_avg = Matrix(grouped_avg[:, [:q10, :q25, :q50, :q75, :q90]])
fanchart!(plt, grouped_avg.bucket, y_avg,
          label="Average Load Error",color=:orange)

# Annual error
grouped_annual = combine(groupby(df, :bucket)) do sub
    q = quantile(sub.annual_error, [0.1, 0.25, 0.5, 0.75, 0.9])
    (; bucket=sub.bucket[1], q10=q[1], q25=q[2], q50=q[3], q75=q[4], q90=q[5])
end
y_annual = Matrix(grouped_annual[:, [:q10, :q25, :q50, :q75, :q90]])
fanchart!(plt, grouped_annual.bucket, y_annual,
          label="Annual Energy Error",color=:green)

display(plt)

###########################################################################################################
###########################################################################################################
#Case study 2: Level of aggregation
###########################################################################################################
###########################################################################################################

#Median of peak errors of MC, different N_L (figure 4 and table 2 in the paper)

###########################################################################################################



plt = plot(xlabel="Number of Buckets B", ylabel="Error")

for (i, N) in enumerate(consumer_sizes)
    df = results_mc[N]

    # Peak error (Upper Peak)
    grouped_peak = combine(groupby(df, :bucket)) do sub
        q = quantile(sub.peak_error, [0.1, 0.25, 0.5, 0.75, 0.9])
        (; bucket=sub.bucket[1], q50=q[3])  # Median only
    end    
    plot!(plt, grouped_peak.bucket, grouped_peak.q50, label="Upper Peak Error (N_L=$N)", lw=2, color=peak_colors[N])

    # Lower Peak error
    grouped_min = combine(groupby(df, :bucket)) do sub
        q = quantile(sub.min_error, [0.1, 0.25, 0.5, 0.75, 0.9])
        (; bucket=sub.bucket[1], q50=q[3])
    end    
    plot!(plt, grouped_min.bucket, grouped_min.q50, label="Lower Peak Error (N_L=$N)", lw=2, color=min_colors[N])
end

# Add threshold lines to the plot
hline!(plt, [0.01], label="Upper Peak Threshold (0.01)", color=:black, linestyle=:dash, lw=2)
#hline!(plt, [0.035], label="Lower Peak Threshold (0.15)", color=:gray, linestyle=:dash, lw=2)
display(plt)

# === Generate Table for Bucket Thresholds ===
# Define error thresholds
thresholds = Dict("Upper Peak" => 0.01)

# Prepare final table
summary_table = DataFrame(N_L=Int[], Upper_Peak=Int[])

for N in consumer_sizes
    df = results_mc[N]

    # Group and compute medians
    grouped = combine(groupby(df, :bucket)) do sub
        (
            bucket=sub.bucket[1],
            upper_median=median(sub.peak_error)
        )
    end

    # Determine bucket counts meeting thresholds
    idx_upper = findfirst(row -> row.upper_median <= thresholds["Upper Peak"], eachrow(grouped))

    min_bucket_upper = isnothing(idx_upper) ? missing : grouped.bucket[idx_upper]
    push!(summary_table, (N, min_bucket_upper))
end

rename!(summary_table, [:N_L, :Upper_Peak])
display(summary_table)

###########################################################################################################
###########################################################################################################
#Case study 3: Application to 76 bus DN
###########################################################################################################
###########################################################################################################

#Post power flow voltage CDF (figure 5 in the paper)

###########################################################################################################
function POLA_network_fluvius_consumers()
    dictionnary_POLA_network = build_mathematical_model(dir, config_file_name; pd = 0.0, qd = 0.0, scale_factor = 1.0)


    #next step would be to build load matrix based upon matching the annual consumption...
    # ... of the load buses of the POLA network with specific fluvius consumers of type 1

    # Step 1: Load the matching CSV file
    matching_file = "matching_POLAnetwork_fluviusconsumers.csv"
    matching_df = CSV.read(matching_file, DataFrame)

    # Extract mappings
    bus_ids_with_load = matching_df[!, "load bus POLA"] .- 1  # correct for offset
    consumer_ids = matching_df[!, "consumer (type 1)"]

    # Step 2: Load total number of buses in the network
    # Assuming you already built the network dictionary
    # (this dictionary was created via your build_mathematical_model function)
    num_buses_total = length(keys(dictionnary_POLA_network["bus"]))

    # Step 3: Define function to load and convert consumer profiles to MW
    function load_consumer_profile_mw(filename::String)
        df = CSV.read(filename, DataFrame, header=false)
        nrows, ncols = size(df)
        data_matrix = Array{Float64}(undef, nrows, ncols)
        for i in 1:nrows, j in 1:ncols
            data_matrix[i, j] = parse(Float64, string(df[i, j]))
        end
        @assert nrows == 96 "Expected 96 rows per day in the consumer file."
        hourly_profile = Float64[]
        for day in 1:ncols
            day_values = data_matrix[:, day]
            day_matrix = reshape(day_values, 4, 24)
            hourly_values = sum(day_matrix, dims=1)
            append!(hourly_profile, vec(hourly_values))
        end
        return hourly_profile ./ 1000  # Convert from kW to MW
    end

    # Step 4: Initialize full load matrix with zeros (for all buses)
    active_load_matrix = zeros(8760, num_buses_total)

    # Directory where consumer CSVs are stored
    load_dir = raw"C:\Users\milan\OneDrive - KU Leuven\2nd master\Thesis\2nd semester\Type1_consumer_Fluvius"

    # Step 5: Fill in the matrix only for buses with matched consumers
    for (consumer_id, bus_id) in zip(consumer_ids, bus_ids_with_load)
        filename = joinpath(load_dir, "Type1consumerLoad$(consumer_id).csv")
        println("Loading consumer $consumer_id → bus $bus_id from $filename")
        active_load_matrix[:, bus_id+1] = load_consumer_profile_mw(filename)
    end

    # Optional: DataFrame with column names as bus IDs
    col_names = ["bus_" * string(i) for i in 1:num_buses_total]
    active_load_df = DataFrame(active_load_matrix, Symbol.(col_names))


    reactive_load_matrix = active_load_matrix .* 0.2  # Assuming a fixed power factor of 0.8 for reactive load
    return dictionnary_POLA_network, active_load_matrix, reactive_load_matrix
end
function convert_to_single_phase_powermodels(pmd_dict::Dict{String,Any})
    pm_dict = Dict{String, Any}()

    pm_dict["baseMVA"] = pmd_dict["baseMVA"]
    pm_dict["per_unit"] = true

    pm_dict["bus"] = Dict{String, Any}()
    pm_dict["branch"] = Dict{String, Any}()
    pm_dict["gen"] = Dict{String, Any}()
    pm_dict["load"] = Dict{String, Any}()
    pm_dict["shunt"] = Dict{String, Any}()
    pm_dict["dcline"] = Dict{String, Any}()
    pm_dict["storage"] = Dict{String, Any}()
    pm_dict["switch"] = Dict{String, Any}()  


    # Buses
    for (bus_id, bus_data) in pmd_dict["bus"]
        pm_dict["bus"][bus_id] = Dict(
            "bus_i" => parse(Int, bus_id),
            "index" => parse(Int, bus_id), 
            "bus_type" => bus_data["bus_type"], 
            "base_kv" => bus_data["vbase"],
            "vm" => haskey(bus_data, "vm") ? bus_data["vm"][1] : 1.0,
            "va" => haskey(bus_data, "va") ? bus_data["va"][1] : 0.0,
            "vmin" => 0.9,
            "vmax" => 1.1
        )
    end
    

    function safe_scalar(x, default)
        if x isa AbstractVector || x isa AbstractMatrix
            return x[1]
        elseif x isa Number
            return x
        else
            return default
        end
    end

    # Branches
    for (br_id, br_data) in pmd_dict["branch"]
        pm_dict["branch"][br_id] = Dict(
            "index" => parse(Int, br_id),  
            "f_bus" => br_data["f_bus"],
            "t_bus" => br_data["t_bus"],
            "br_r" => br_data["br_r"][1,1],
            "br_x" => br_data["br_x"][1,1],
            "br_status" => br_data["br_status"],
            "angmin" => -2π,
            "angmax" => 2π,
            "tap"   => safe_scalar(get(br_data, "tap", 1.0), 1.0),
            "shift" => safe_scalar(get(br_data, "shift", 0.0), 0.0),
            "g_fr"  => safe_scalar(get(br_data, "g_fr", [0.0]), 0.0),
            "b_fr"  => safe_scalar(get(br_data, "b_fr", [0.0]), 0.0),
            "g_to"  => safe_scalar(get(br_data, "g_to", [0.0]), 0.0),
            "b_to"  => safe_scalar(get(br_data, "b_to", [0.0]), 0.0),
            "tr"    => safe_scalar(get(br_data, "tr", [1.0]), 1.0),
            "ti"    => safe_scalar(get(br_data, "ti", [0.0]), 0.0),
            "tm"    => safe_scalar(get(br_data, "tm", [1.0]), 1.0)
        )
    end

    

    # Loads
    for (ld_id, ld_data) in pmd_dict["load"]
        pm_dict["load"][ld_id] = Dict(
            "load_bus" => ld_data["load_bus"],
            "status" => ld_data["status"],
            "pd" => ld_data["pd"][1],  
            "qd" => ld_data["qd"][1]
        )
    end

    # Generator
    if haskey(pmd_dict, "gen")
        for (gen_id, gen_data) in pmd_dict["gen"]
            gen_status = haskey(gen_data, "gen_status") ? gen_data["gen_status"] :
                        haskey(gen_data, "status") ? gen_data["status"] : 1

            pm_dict["gen"][gen_id] = Dict(
                "index" => parse(Int, gen_id),
                "gen_bus" => gen_data["gen_bus"],
                "pg" => gen_data["pg"][1],
                "qg" => gen_data["qg"][1],
                "pmax" => gen_data["pmax"][1],
                "pmin" => gen_data["pmin"][1],
                "qmax" => gen_data["qmax"][1],
                "qmin" => gen_data["qmin"][1],
                "gen_status" => gen_status,  # <-- DIT is wat PowerModels verwacht!
                "vg" => 1.0,
                "cost" => [0.0, 0.0],
                "ncost" => 2
            )
        end
    end




    return pm_dict
end

POLA_dict, active_load_matrix, reactive_load_matrix = POLA_network_fluvius_consumers()
POLA_dict= convert_to_single_phase_powermodels(POLA_dict)


n_scenarios = size(active_load_matrix, 1)
n_buses = length(POLA_dict["bus"])
vm_matrix = zeros(n_scenarios, n_buses)
termination_statuses = Vector{String}(undef, n_scenarios)


# Loop over all load scenarios
for scenario_idx in 1:n_scenarios
    # Create a fresh copy of the network for this scenario
    println("Running powerflow for scenario $scenario_idx")
    network = deepcopy(POLA_dict)

    # Keep track of which buses already had their load updated
    updated_buses = Set{Int}()

    # Update active and reactive loads according to the scenario
    for (load_id, load_data) in network["load"]
        bus_idx = load_data["load_bus"]  

        if !(bus_idx in updated_buses)
            # If the bus was not yet updated, fill with scenario load
            load_data["pd"] = active_load_matrix[scenario_idx, bus_idx]
            load_data["qd"] = reactive_load_matrix[scenario_idx, bus_idx]

            push!(updated_buses, bus_idx)  # Mark this bus as updated
        else
            # If the bus was already updated, set pd and qd to zero
            load_data["pd"] = 0.0
            load_data["qd"] = 0.0
        end
    end


    # Solve the power flow
    result = PowerModels.solve_pf(network, ACPPowerModel, Ipopt.Optimizer)
    
    # Save the termination status (always)
    termination_statuses[scenario_idx] = string(result["termination_status"])
    # Check if the power flow was successful
    if result["termination_status"] == PowerModels.LOCALLY_SOLVED
        bus_solution = result["solution"]["bus"]
        for (bus_id, bus_data) in bus_solution
            bus_index = parse(Int, bus_id)
            vm_matrix[scenario_idx, bus_index] = bus_data["vm"]
        end
    else
        println("Power flow failed for scenario $scenario_idx")
    end
end

#save the vm_matrix to a csv file
CSV.write("vm_matrix_postpowerflowevaluation.csv", DataFrame(vm_matrix), header=["bus_" * string(i) for i in 1:n_buses])

# compute the reduced scenarios for the POLA network (with assigned fluvius consumers)
reduced_info = casestudy4_reducedscenarios()

# Initialize an empty matrix for reconstructed voltages
vm_matrix_reconstructed = Array{Float64}(undef, 0, n_buses)  # Start empty

for row in eachrow(reduced_info)
    original_hour = row[Symbol("hour of the year")]
    final_weight = row[Symbol("final_weight")]

    # How many hours this reduced scenario should represent
    num_hours = round(Int, final_weight * 8760)
    
    # Fetch the voltage profile at the original hour
    voltage_row = vm_matrix[original_hour, :]  # 1 row, all buses
    
    # Repeat the row num_hours times and append to reconstructed matrix
    repeated_rows = repeat(reshape(voltage_row, 1, :), num_hours, 1)
    vm_matrix_reconstructed = vcat(vm_matrix_reconstructed, repeated_rows)
end

# Check final size
println("Final reconstructed vm matrix size: ", size(vm_matrix_reconstructed))
current_size = size(vm_matrix_reconstructed, 1)
missing_rows = 8760 - current_size

if missing_rows > 0
    # Add missing rows by repeating the last row
    last_row = vm_matrix_reconstructed[end, :]
    repeated_rows = repeat(reshape(last_row, 1, :), missing_rows, 1)
    vm_matrix_reconstructed = vcat(vm_matrix_reconstructed, repeated_rows)
elseif missing_rows < 0
    # Remove extra rows
    vm_matrix_reconstructed = vm_matrix_reconstructed[1:end+missing_rows, :]
end
############################
#BUILD CDFs to check the post powerflow evaluation of the reduced scenarios
############################

# Flatten and sort the voltage matrices
vm_values_original = vec(vm_matrix)
vm_values_reduced = vec(vm_matrix_reconstructed)

# Number of probability points
n = length(vm_values_original)  # could also set n = 1000 for smoother plot
p_vals = range(0, 1; length=n)

# Compute quantiles (inverse CDFs)
quant_orig = quantile(vm_values_original, p_vals)
quant_red = quantile(vm_values_reduced, p_vals)


# Plot the original and reduced CDF curves
plot(quant_orig, p_vals,
     label="All Load Scenarios", color=:blue, lw=2)
plot!(quant_red, p_vals,
      label="Reduced Scenarios", color=:orange, lw=2)


# Add labels and title
ylabel!("Cumulative Probability")
xlabel!(L"Voltage Magnitude $v_m$ (p.u.)")