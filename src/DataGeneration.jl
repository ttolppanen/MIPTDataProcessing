export generate_data

function get_alg_types()
    return [:exact, :krylov, :mps]
end

function generate_data(; kwargs...)
    if kwargs[:alg] == :exact
        calc_exact(; kwargs[])
    else if kwargs[:alg] == :krylov
        calc_krylov(; d = kwargs[:d], L = kwargs[:L], dt = kwargs[:dt], t = kwargs[:t], state = kwargs[:state], k = kwargs[:k], p = kwargs[:p], traj = kwargs[:traj])
    else if kwargs[:alg] == :mps
    else
        println("error in alg")
        # throw error
    end
end

function calc_exact()
end

function calc_krylov(; d, L, dt, t, state, k, p, traj)
    state = eval(Meta.parse(state))
    H = bosehubbard(d, L)
    effect! = measuresitesrandomly!(L, nop(d), p)
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