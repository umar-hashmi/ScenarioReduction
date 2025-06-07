# This code handles the scenario generation for the load profile, a test
# Pkg.add("Distributions")
using Distributions
using LinearAlgebra
using Random
using PowerModels

# Function to generate scenarios
function generate_scenarios(mu::Vector{Float64}, sigma::Vector{Float64}, corr_matrix::Matrix{Float64}, J::Int, epsilon::Float64=1e-10)
    # Number of intervals
    S = length(mu)
    # Create the covariance matrix Σ
    Σ = diagm(sigma) * corr_matrix * diagm(sigma)

    # Ensure Σ is symmetric by averaging it with its transpose
    Σ = (Σ + Σ') / 2

    # Add a small error to the covariance matrix for numerical stability
    Σ += epsilon * I(S)

    # Perform Cholesky decomposition
    L = cholesky(Σ).L

    # Generate random numbers from a standard normal distribution
    rnd = randn(S, J)

    # Generate scenarios using the formula X = μ + LZ
    scenarios = (mu .+ L * rnd)'

    return scenarios
end

# Function to generate load profiles (scenarios)
function scen_gen(load_profiles, S, N, forecast_err)
    mu = collect(load_profiles)  # Mean vector from load profiles
    sigma = mu * forecast_err  # Standard deviation as fraction of the mean
    corr_matrix = fill(forecast_err, N, N) + I(N)*(1 - forecast_err)  # Simple correlation matrix
    load_sc = generate_scenarios(mu, sigma, corr_matrix, S)
    return load_sc
end

# Main function to generate load scenarios
#for IEEE33 kW to MW conversion (not correct in .m file)
function load_scenario_generation_IEEE33(data::Dict{}, num_scenarios::Int)
    
    # Extract data
    bus_data = data["bus"]  # All buses as a dictionary
    bus_loads = get(data, "load", Dict())  # Load data, default to empty Dict if not present

    N = length(bus_data)  # Number of buses

    # Step 1: Correct mapping of bus IDs to array indices
    bus_ids = sort(parse.(Int, collect(keys(bus_data))))  # Get sorted bus IDs
    bus_index_map = Dict(string(bus_id) => i for (i, bus_id) in enumerate(bus_ids))  # Map bus ID to index

    # Step 2: Initialize load vectors with zeros
    mu_Pd = zeros(N)
    mu_Qd = zeros(N)

    # Step 3: Populate the vectors with load data where available
    for (key, load) in bus_loads
        bus_id = string(load["load_bus"])  # Use the "load_bus" to find the correct bus
        if haskey(bus_index_map, bus_id)  # Check if bus ID exists in the mapping
            idx = bus_index_map[bus_id]  # Get the correct index
            mu_Pd[idx] = load["pd"] / 1e3  # Active load in MW
            mu_Qd[idx] = load["qd"] / 1e3  # Reactive load in MW
        end
    end

    # Generate load scenarios
    forecast_err = 0.3
    active_load_scenarios = scen_gen(mu_Pd, num_scenarios, N, forecast_err)
    reactive_load_scenarios = scen_gen(mu_Qd, num_scenarios, N, forecast_err)

    return active_load_scenarios, reactive_load_scenarios
end
#for ENWL and other networks where load is just in MW (correct for pu)
function load_scenario_generation(data::Dict{}, num_scenarios::Int)
    
    # Extract data
    bus_data = data["bus"]  # All buses as a dictionary
    bus_loads = get(data, "load", Dict())  # Load data, default to empty Dict if not present

    N = length(bus_data)  # Number of buses

    # Step 1: Correct mapping of bus IDs to array indices
    bus_ids = sort(parse.(Int, collect(keys(bus_data))))  # Get sorted bus IDs
    bus_index_map = Dict(string(bus_id) => i for (i, bus_id) in enumerate(bus_ids))  # Map bus ID to index

    # Step 2: Initialize load vectors with zeros
    mu_Pd = zeros(N)
    mu_Qd = zeros(N)

    # Step 3: Populate the vectors with load data where available
    for (key, load) in bus_loads
        bus_id = string(load["load_bus"])  # Use the "load_bus" to find the correct bus
        if haskey(bus_index_map, bus_id)  # Check if bus ID exists in the mapping
            idx = bus_index_map[bus_id]  # Get the correct index
            mu_Pd[idx] = load["pd"]   # Active load in MW
            mu_Qd[idx] = load["qd"]   # Reactive load in MW
        end
    end

    # Generate load scenarios
    forecast_err = 0.3  # 10% forecast error
    active_load_scenarios = scen_gen(mu_Pd, num_scenarios, N, forecast_err)
    reactive_load_scenarios = scen_gen(mu_Qd, num_scenarios, N, forecast_err)

    return active_load_scenarios, reactive_load_scenarios
end