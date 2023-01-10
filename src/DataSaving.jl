# using HDF5

function saveh5(filename, obsrv_data, msr_prob; observable, simulation_param...)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    h5open(path_to_file, "cw") do file
        for g_mipt in file
            if same_param(g_mipt, simulation_param)
                if haskey(g_mipt, "p = $msr_prob")
                    println("Data already exists for p = $msr_prob")
                    return
                else
                    g_p = new_p_group(g_mipt, msr_prob)
                    add_data_to_p_group(g_p, obsrv_data; observable)
                end
                return
            end
        end
        g_mipt = new_mipt_group(file, simulation_param...)
        g_p = new_p_group(g_mipt, msr_prob)
        add_data_to_p_group(g_p, obsrv_data; observable)
    end
end

function new_mipt_group(file; simulation_param...)
    group_name = "mipt_" * string(length(file) + 1)
    g_mipt = create_group(file, group_name)
    for (key, value) in simulation_param
        attributes(g_mipt)["$key"] = value
    end
    return g_mipt
end

function new_p_group(g_mipt, msr_prob)
    g_p = create_group(g_mipt, "p = $(msr_prob)")
    attributes(g_p)["p"] = msr_prob   
    return g_p
end

function add_data_to_p_group(g_p, obsrv_data; observable)
    data_size = size(obsrv_data) # 1 timesteps, 2 trajectories
    data = create_dataset(g_p, string(observable), Float64, ((data_size[1], data_size[2]), (data_size[1], -1)), chunk = (data_size[1], 1))
    data[:, :] = obsrv_data
    attributes(data)["end_mean"] = mean(obsrv_data[end, :])
end

