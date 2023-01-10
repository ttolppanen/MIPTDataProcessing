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
    state0 = eval(Meta.parse(get_initial_states()[state][:exact]))
    H = bosehubbard(d, L; w, U, J)
    U = exp(-im * dt * Matrix(H))
    effect! = eval(Meta.parse(get_measurements()[msr][:exact]))
    r_f() = exactevolve(state0, U, dt, t; effect!)
    r = solvetrajectories(r_f, traj)
    for observable in observables
        obsrv_data = zeros(length(r[1]), length(r))
        for traj_i in eachindex(r)
            for timestep_i in eachindex(r[traj_i])
                s = r[traj_i][timestep_i] # this gets used in get_observable
                obsrv_data[timestep_i, traj_i] = eval(Meta.parse(get_observables()[observable][:exact]))
            end
        end
        saveh5(filename, obsrv_data, p; d, L, dt, t, state, msr, w, U, J, observable, alg)
    end
end

function calc_krylov(filename; d, L, dt, t, state, k, p, msr, traj, w, U, J, observables, alg)
    println("JOO")
    state0 = eval(Meta.parse(get_initial_states()[state][:exact]))
    H = bosehubbard(d, L; w, U, J)
    effect! = eval(Meta.parse(get_measurements()[msr][:exact]))
    r_f() = krylovevolve(state0, H, dt, t, k; effect!)
    r = solvetrajectories(r_f, traj)
    for observable in observables
        obsrv_data = zeros(length(r[1]), length(r))
        for traj_i in eachindex(r)
            for timestep_i in eachindex(r[traj_i])
                s = r[traj_i][timestep_i] # this gets used in get_observable
                obsrv_data[timestep_i, traj_i] = eval(Meta.parse(get_observables()[observable][:exact]))
            end
        end
        saveh5(filename, obsrv_data, p; d, L, dt, t, state, k, msr, w, U, J, observable, alg)
    end
end

function calc_mps(filename; d, L, dt, t, state, trotter_order, p, msr, traj, w, U, J, observables, alg, ITensors_apply_kwargs...)
    state0 = eval(Meta.parse(get_initial_states()[state][:mps]))
    gates = bosehubbardgates(siteinds(state0), dt; k=trotter_order)
    effect! = eval(Meta.parse(get_measurements()[msr][:mps]))
    r_f() = mpsevolve(state0, gates, dt, t; effect!, ITensors_apply_kwargs...)
    r = solvetrajectories(r_f, traj)
    for observable in observables
        obsrv_data = zeros(length(r[1]), length(r))
        for traj_i in eachindex(r)
            for timestep_i in eachindex(r[traj_i])
                s = r[traj_i][timestep_i] # this gets used in get_observable
                obsrv_data[timestep_i, traj_i] = eval(Meta.parse(get_observables()[observable][:mps]))
            end
        end
        saveh5(filename, obsrv_data, p; d, L, dt, t, state, trotter_order, msr, w, U, J, observable, alg)
    end
end