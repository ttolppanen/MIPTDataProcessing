# using QuantumStates
# using QuantumOperators
# using QuantumTimeEvolution

export generate_data

function get_msr(msr)
    measurement = Dict(
        "std_n" => (;exact = "measuresitesrandomly!(L, nop(d), p)", mps = "asdasd")
    )
    return measurement[msr]
end

function generate_data(; kwargs...)
    if kwargs[:alg] == :exact
        calc_exact(; kwargs[])
    else if kwargs[:alg] == :krylov
        calc_krylov(; d = kwargs[:d], L = kwargs[:L], dt = kwargs[:dt], t = kwargs[:t], state = kwargs[:state], k = kwargs[:k], p = kwargs[:p], msr_type = kwargs[:msr_type], traj = kwargs[:traj])
    else if kwargs[:alg] == :mps
    else
        println("error in alg")
        # throw error
    end
end

function calc_exact()
end

function calc_krylov(; d, L, dt, t, state, k, p, msr_type, traj)
    state = eval(Meta.parse(state))
    H = bosehubbard(d, L)
    effect! = eval(Meta.parse(get_msr(msr_type)[:exact]))
    r_f() = krylovevolve(state, H, dt, t, k; effect!)
    r = solvetrajectories(r_f, traj)
    ent_data = zeros(length(r[1]), length(r))
    for traj_i in eachindex(r) 
        for timestep_i in eachindex(r[traj_i])
            s = r[traj_i][timestep_i]
            ent_data[timestep_i, traj_i] = entanglement(d, L, s, Int(floor(L/2)))
        end
    end
    saveh5(ent_data, p; d, L, dt, t, k, initial_state = "zeroone(d, L)", w=1, U=1, J=1, msr_type = "std_n")
end