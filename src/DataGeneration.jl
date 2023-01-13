# using QuantumStates
# using QuantumOperators
# using QuantumTimeEvolution

export generate_data

function generate_data(filename; kwargs...)
    if kwargs[:alg] == :exact
        calc_exact(filename; kwargs...)
    elseif kwargs[:alg] == :krylov 
        calc_krylov(filename; kwargs...)
    elseif kwargs[:alg] == :mps
        calc_mps(filename; kwargs...)
    else
        println("error in alg")
        # throw error
    end
end

function calc_exact(filename; d, L, dt, t, state, probabilities, measurement, traj, w, U, J, observables, alg)
    sp = simulation_parameters_to_kwargs(; d, L, dt, t, state, measurement, w, U, J)
    state0 = get_initial_states()[state][:exact](; sp...)
    H = bosehubbard(d, L; w, U, J)
    expU = exp(-im * dt * Matrix(H))
    meffect! = get_measurements()[measurement][:exact](; sp...)
    for p in probabilities
        effect!(state) = meffect!(state, p)
        r_f() = exactevolve(state0, expU, dt, t; effect!)
        r = solvetrajectories(r_f, traj)
        for observable in observables
            #get correct num of traj
            #if c_t < 1 skip
            #r = solvetrajectories(r_f, traj_real)
            obsrv_data = zeros(length(r[1]), length(r))
            for traj_i in eachindex(r)
                for timestep_i in eachindex(r[traj_i])
                    s = r[traj_i][timestep_i] # this gets used in get_observable
                    obsrv_data[timestep_i, traj_i] = get_observables()[observable][:exact](s; sp...)
                end
            end
            saveh5(filename, obsrv_data, p; d, L, dt, t, state, measurement, w, U, J, observable, alg)
        end
    end
end

function calc_krylov(filename; d, L, dt, t, state, k, probabilities, measurement, traj, w, U, J, observables, alg)
    sp = simulation_parameters_to_kwargs(; d, L, dt, t, state, k, measurement, w, U, J)
    state0 = get_initial_states()[state][:exact](; sp...)
    H = bosehubbard(d, L; w, U, J)
    meffect! = get_measurements()[measurement][:exact](; sp...)
    for p in probabilities
        effect!(state) = meffect!(state, p)
        r_f() = krylovevolve(state0, H, dt, t, k; effect!)
        r = solvetrajectories(r_f, traj)
        for observable in observables
            obsrv_data = zeros(length(r[1]), length(r))
            for traj_i in eachindex(r)
                for timestep_i in eachindex(r[traj_i])
                    s = r[traj_i][timestep_i] # this gets used in get_observable
                    obsrv_data[timestep_i, traj_i] = get_observables()[observable][:exact](s; sp...)
                end
            end
            saveh5(filename, obsrv_data, p; d, L, dt, t, state, k, measurement, w, U, J, observable, alg)
        end
    end
end

function calc_mps(filename; d, L, dt, t, state, trotter_order, probabilities, measurement, traj, w, U, J, observables, alg, ITensors_apply_kwargs...)
    sp = simulation_parameters_to_kwargs(; d, L, dt, t, state, trotter_order, measurement, w, U, J)
    state0 = get_initial_states()[state][:mps](; sp...)
    gates = bosehubbardgates(siteinds(state0), dt; k=trotter_order)
    meffect! = get_measurements()[measurement][:mps](ITensors_apply_kwargs; state0, sp...)
    for p in probabilities
        effect!(state) = meffect!(state, p)
        r_f() = mpsevolve(state0, gates, dt, t; effect!, ITensors_apply_kwargs...)
        r = solvetrajectories(r_f, traj)
        for observable in observables
            obsrv_data = zeros(length(r[1]), length(r))
            for traj_i in eachindex(r)
                for timestep_i in eachindex(r[traj_i])
                    s = r[traj_i][timestep_i] # this gets used in get_observable
                    obsrv_data[timestep_i, traj_i] = get_observables()[observable][:mps](s; sp...)
                end
            end
            saveh5(filename, obsrv_data, p; d, L, dt, t, state, trotter_order, measurement, w, U, J, observable, alg, ITensors_apply_kwargs...)
        end
    end
end

function simulation_parameters_to_kwargs(; kwargs...)
    return kwargs
end