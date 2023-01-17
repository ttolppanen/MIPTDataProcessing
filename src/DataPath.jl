export setdatapath
export getdatapath
export get_num_of_traj
export get_min_existing_traj

function setdatapath(path)
    open(joinpath(@__DIR__, "..", "datapath.txt"), "w") do file
        write(file, path)
    end
end

function getdatapath()
    open(joinpath(@__DIR__, "..", "datapath.txt"), "r") do file
        read(file, String)
    end
end

function get_num_of_traj(filename, msr_prob, observable; simulation_param...)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    try
        h5open(path_to_file, "r") do file
            for g_mipt in file
                if same_param(g_mipt, simulation_param)
                    return size(g_mipt["p = $msr_prob"][string(observable)])[2] # 1 timesteps, 2 trajectories
                end
            end
            return 0;
        end
    catch
        return 0;
    end
end

function get_min_existing_traj(filename, msr_prob, observables; simulation_param...)
    out = 0
    for i in eachindex(observables)
        observable = observables[i]
        if i == 1
            out = get_num_of_traj(filename, msr_prob, observable; simulation_param...)
        else
            traj = get_num_of_traj(filename, msr_prob, observable; simulation_param...)
            if (traj > out)
                out = traj 
            end
        end
    end
    return out
end
