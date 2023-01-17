# using HDF5
# using Statistics

function saveh5(filename, obsrv_data, msr_prob, observable; simulation_param...)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    h5open(path_to_file, "cw") do file
        for g_mipt in file
            if same_param(g_mipt, simulation_param)
                if haskey(g_mipt, "p = $msr_prob")
                    g_p = g_mipt["p = $msr_prob"]
                    if haskey(g_p, string(observable))
                        dset_obsrv = g_p[string(observable)]
                        add_data_to_obsrv_group(dset_obsrv, obsrv_data)
                    else
                        create_data_to_p_group(g_p, obsrv_data, observable)
                    end
                    return
                else
                    g_p = new_p_group(g_mipt, msr_prob)
                    create_data_to_p_group(g_p, obsrv_data, observable)
                end
                return
            end
        end
        g_mipt = new_mipt_group(file; simulation_param...)
        g_p = new_p_group(g_mipt, msr_prob)
        create_data_to_p_group(g_p, obsrv_data, observable)
    end
end

function new_mipt_group(file; simulation_param...)
    group_name = "mipt_" * string(length(file) + 1)
    g_mipt = create_group(file, group_name)
    for (key, value) in simulation_param
        attributes(g_mipt)["$key"] = isa(value, Symbol) ? string(value) : value
    end
    return g_mipt
end

function new_p_group(g_mipt, msr_prob)
    g_p = create_group(g_mipt, "p = $(msr_prob)")
    attributes(g_p)["p"] = msr_prob   
    return g_p
end

function create_data_to_p_group(g_p, obsrv_data, observable)
    data_size = size(obsrv_data) # 1 timesteps, 2 trajectories
    data = create_dataset(g_p, string(observable), Float64, ((data_size[1], data_size[2]), (data_size[1], -1)), chunk = (data_size[1], 1))
    data[:, :] = obsrv_data
    attributes(data)["end_mean"] = mean(obsrv_data[end, :])
    attributes(data)["end_standard_deviation"] = std(obsrv_data[end, :])
end

function add_data_to_obsrv_group(dset_obsrv, obsrv_data)
    dsetsize = size(dset_obsrv) # 1 timesteps, 2 trajectories
    timesteps = dsetsize[1]
    traj = dsetsize[2]
    new_traj = size(obsrv_data)[2]
    HDF5.set_extent_dims(dset_obsrv, (timesteps, traj + new_traj))
    dset_obsrv[:, (traj + 1):(traj + new_traj)] = obsrv_data
    end_mean = mean(dset_obsrv[end, :])
    end_std = std(dset_obsrv[end, :])
    delete_attribute(dset_obsrv, "end_mean")
    write_attribute(dset_obsrv, "end_mean", end_mean)
    delete_attribute(dset_obsrv, "end_standard_deviation")
    write_attribute(dset_obsrv, "end_standard_deviation", end_std)
end

function same_param(g, params)
    if issetequal(keys(attributes(g)), String.(keys(params)))
        attribute_values = [read_attribute(g, key) for key in keys(attributes(g))]
        param_values = [isa(params[Symbol(key)], Symbol) ? string(params[Symbol(key)]) : params[Symbol(key)] for key in keys(attributes(g))]
        if issetequal(attribute_values, param_values)
            return true
        end
    end
    return false
end