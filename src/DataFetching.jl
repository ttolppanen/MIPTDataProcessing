# using HDF5

export get_probabilities
export get_data_end_mean
export get_attributes_string

function get_probabilities(filename, groupname)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    h5open(path_to_file, "r") do file
        g_mipt = file[groupname]
        return get_probabilities(g_mipt)
    end
end
function get_probabilities(group::HDF5.Group)
    out = []
    for g_p in group
        push!(out, read_attribute(g_p, "p"))
    end
    return out
end

function get_data_end_mean(filename, groupname, probabilities, observer::Symbol)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    h5open(path_to_file, "r") do file
        g_mipt = file[groupname]
        return get_data_end_mean(g_mipt, probabilities, observer)
    end
end
function get_data_end_mean(group::HDF5.Group, probabilities, observer::Symbol)
    out = []
    for p in probabilities
        e_mean = read_attribute(group["p = $p"][string(observer)], "end_mean")
        push!(out, e_mean)
    end
    return out
end

function get_attributes_string(filename, groupname)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    h5open(path_to_file, "r") do file
        g_mipt = file[groupname]
        return get_attributes_string(g_mipt)
    end
end
function get_attributes_string(group::HDF5.Group)
    out = ""
    for attr_key in keys(attributes(group))
        value = string(read_attribute(group, attr_key))
        out *= "$attr_key = $value, "
    end
    return out[1:end-2]
end