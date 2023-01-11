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

function calc_exact(filename; d, L, dt, t, state, p, msr, traj, w, U, J, observables, alg)
    sp = simulation_parameters_to_kwargs(; d, L, dt, t, state, p, msr, w, U, J)
    state0 = get_initial_states()[state][:exact](; sp...)
    H = bosehubbard(d, L; w, U, J)
    expU = exp(-im * dt * Matrix(H))
    effect! = get_measurements()[msr][:exact](; sp...)
    r_f() = exactevolve(state0, expU, dt, t; effect!)
    r = solvetrajectories(r_f, traj)
    for observable in observables
        obsrv_data = zeros(length(r[1]), length(r))
        for traj_i in eachindex(r)
            for timestep_i in eachindex(r[traj_i])
                s = r[traj_i][timestep_i] # this gets used in get_observable
                obsrv_data[timestep_i, traj_i] = get_observables()[observable][:exact](s; sp...)
            end
        end
        saveh5(filename, obsrv_data, p; d, L, dt, t, state, msr, w, U, J, observable, alg)
    end
end

function calc_krylov(filename; d, L, dt, t, state, k, p, msr, traj, w, U, J, observables, alg)
    sp = simulation_parameters_to_kwargs(; d, L, dt, t, state, k, p, msr, w, U, J)
    state0 = get_initial_states()[state][:exact](; sp...)
    H = bosehubbard(d, L; w, U, J)
    effect! = get_measurements()[msr][:exact](; sp...)
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
        saveh5(filename, obsrv_data, p; d, L, dt, t, state, k, msr, w, U, J, observable, alg)
    end
end

function calc_mps(filename; d, L, dt, t, state, trotter_order, p, msr, traj, w, U, J, observables, alg, ITensors_apply_kwargs...)
    sp = simulation_parameters_to_kwargs(; d, L, dt, t, state, trotter_order, p, msr, w, U, J)
    state0 = get_initial_states()[state][:mps](; sp...)
    gates = bosehubbardgates(siteinds(state0), dt; k=trotter_order)
    effect! = get_measurements()[msr][:mps](ITensors_apply_kwargs; state0, sp...)
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
        saveh5(filename, obsrv_data, p; d, L, dt, t, state, trotter_order, msr, w, U, J, observable, alg, ITensors_apply_kwargs...)
    end
end

function simulation_parameters_to_kwargs(; kwargs...)
    return kwargs
end