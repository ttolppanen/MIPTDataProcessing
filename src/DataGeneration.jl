# using QuantumStates
# using QuantumOperators
# using QuantumTimeEvolution

export generate_data

function generate_data(filename; kwargs...)
    check_version_number(filename)
    if kwargs[:alg] == :exact
        calc_exact(filename; kwargs...)
    elseif kwargs[:alg] == :krylov 
        calc_krylov(filename; kwargs...)
    elseif kwargs[:alg] == :mps
        calc_mps(filename; kwargs...)
    else
        throw(error("No alg provided. Possible alg are :exact, :krylov, :mps"))
    end
    add_version_number(filename)
end

function calc_exact(filename; d, L, dt, t, state, probabilities, measurement, traj, w, U, J, observables, alg, save_before_effect)
    sp = simulation_parameters_to_kwargs(; d, L, dt, t, state, measurement, w, U, J, alg, save_before_effect)
    state0 = get_initial_states()[state][:exact](; sp...)
    H = bosehubbard(d, L; w, U, J)
    expU = exp(-im * dt * Matrix(H))
    meffect! = get_measurements()[measurement][:exact](; sp...)
    for p in probabilities
        traj_to_solve = trajectories_for_solving(filename, traj, p, observables; sp...)
        if traj_to_solve <= 0 continue end # skip this loop if the data exists already
        obs_f = [get_observables()[obs][:exact](; sp...) for obs in observables]
        effect!(state) = meffect!(state, p)
        r_f() = exactevolve(state0, expU, dt, t, obs_f...; effect!, save_before_effect)
        r = solvetrajectories(r_f, traj_to_solve; paral = :distributed)
        for (obs_i, observable) in enumerate(observables)
            traj_to_calculate = trajectories_for_calculating_observables(filename, traj, p, observable; sp...)
            if traj_to_calculate <= 0 continue end # skip this loop if the data exists already
            obsrv_data = extract_observable_data(r, obs_i, traj_to_calculate)
            saveh5(filename, obsrv_data, p, observable; sp...)
        end
    end
end

function calc_krylov(filename; d, L, dt, t, state, k, probabilities, measurement, traj, w, U, J, observables, alg, save_before_effect)
    sp = simulation_parameters_to_kwargs(; d, L, dt, t, state, k, measurement, w, U, J, alg, save_before_effect)
    state0 = get_initial_states()[state][:exact](; sp...)
    H = bosehubbard(d, L; w, U, J)
    meffect! = get_measurements()[measurement][:exact](; sp...)
    for p in probabilities
        traj_to_solve = trajectories_for_solving(filename, traj, p, observables; sp...)
        if traj_to_solve <= 0 continue end # skip this loop if the data exists already
        obs_f = [get_observables()[obs][:exact](; sp...) for obs in observables]
        effect!(state) = meffect!(state, p)
        r_f() = krylovevolve(state0, H, dt, t, k, obs_f...; effect!, save_before_effect)
        r = solvetrajectories(r_f, traj_to_solve; paral = :distributed)
        for (obs_i, observable) in enumerate(observables)
            traj_to_calculate = trajectories_for_calculating_observables(filename, traj, p, observable; sp...)
            if traj_to_calculate <= 0 continue end # skip this loop if the data exists already
            obsrv_data = extract_observable_data(r, obs_i, traj_to_calculate)
            saveh5(filename, obsrv_data, p, observable; sp...)
        end
    end
end

function calc_mps(filename; d, L, dt, t, state, trotter_order, probabilities, measurement, traj, w, U, J, observables, alg, save_before_effect, ITensors_apply_kwargs...)
    sp = simulation_parameters_to_kwargs(; d, L, dt, t, state, trotter_order, measurement, w, U, J, alg, save_before_effect)
    state0 = get_initial_states()[state][:mps](; sp...)
    gates = bosehubbardgates(siteinds(state0), dt; k=trotter_order)
    meffect! = get_measurements()[measurement][:mps](ITensors_apply_kwargs; state0, sp...)
    for p in probabilities
        traj_to_solve = trajectories_for_solving(filename, traj, p, observables; sp..., ITensors_apply_kwargs...)
        if traj_to_solve <= 0 continue end # skip this loop if the data exists already
        obs_f = [get_observables()[obs][:mps](; sp...) for obs in observables]
        effect!(state) = meffect!(state, p)
        r_f() = mpsevolve(state0, gates, dt, t, obs_f...; effect!, save_before_effect, ITensors_apply_kwargs...)
        r = solvetrajectories(r_f, traj_to_solve; paral = :distributed)
        for (obs_i, observable) in enumerate(observables)
            traj_to_calculate = trajectories_for_calculating_observables(filename, traj, p, observable; sp..., ITensors_apply_kwargs...)
            if traj_to_calculate <= 0 continue end # skip this loop if the data exists already
            obsrv_data = extract_observable_data(r, obs_i, traj_to_calculate)
            saveh5(filename, obsrv_data, p, observable; sp..., ITensors_apply_kwargs...)
        end
    end
end

function check_version_number(filename)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    if isfile(path_to_file)
        h5open(path_to_file, "cw") do file
            if haskey(attributes(file), "version")
                file_version = read_attribute(file, "version")
                package_version = package_version_number()
                error_text = "Package version is $package_version but file version is $file_version"
                file_version = split(file_version, ".")
                package_version = split(package_version, ".")
                if file_version[1] != package_version[1] || file_version[2] != package_version[2]
                    throw(error(error_text))
                end
            else
                throw(error("File exists but doesn't have a version number"))
            end
        end
    end
end

function add_version_number(filename)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    h5open(path_to_file, "cw") do file
        if (!haskey(attributes(file), "version"))
            attributes(file)["version"] = package_version_number()
        end
    end
end

function simulation_parameters_to_kwargs(; kwargs...)
    return kwargs
end

function trajectories_for_solving(filename, traj, p, observables; sp...)
    traj_to_solve = traj - get_min_existing_traj(filename, p, observables; sp...)
    if traj_to_solve <= 0
        println("Trajectories already exist for all observers in p = " * string(p))
    else
        println("Solving " * string(traj_to_solve) * " trajectories for p = " * string(p))
    end
    return traj_to_solve
end

function trajectories_for_calculating_observables(filename, traj, p, observable; sp...)
    traj_to_calculate = traj - get_num_of_traj(filename, p, observable; sp...)
    if traj_to_calculate <= 0
        println("Trajectories already exist for p = " * string(p) * " observable = " * string(observable))
    end
    return traj_to_calculate
end

function extract_observable_data(result, obs_i, traj_to_calculate)
    num_of_timesteps = length(result[1][obs_i, :])
    out = zeros(num_of_timesteps, traj_to_calculate)
    for traj_i in 1:traj_to_calculate
        for timestep_i in 1:num_of_timesteps
            out[timestep_i, traj_i] = result[traj_i][obs_i, timestep_i]
        end
    end
    return out
end