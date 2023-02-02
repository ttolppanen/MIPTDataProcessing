# using HDF5

export get_groups
export get_groups_with_param
export get_probabilities
export get_data
export get_data_end_mean
export get_attributes_string
export get_attribute

function get_groups(filename)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    h5open(path_to_file, "r") do file
        return keys(file)
    end
end

function get_groups_with_param(filename; params...)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    h5open(path_to_file, "r") do file
        for g in file
            if all([haskey(attributes(g), string(key)) && read_attribute(g, string(key)) == value for (key, value) in params])
                println(HDF5.name(g) * " | " * get_attributes_string(g; addnewline = false))
            end
        end
    end
end

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

function get_data(filename::String, groupname::String, probability::Real, observer::String, trajectories = :all)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    h5open(path_to_file, "r") do file
        g_mipt = file[groupname]
        return get_data(g_mipt, probability, observer, trajectories)
    end
end
function get_data(group::HDF5.Group, probability::Real, observer::String, trajectories = :all)
    if (trajectories == :all)
        return read(group["p = $probability"][observer])
    else
        return group["p = $probability"][observer][1:end, trajectories]
    end
end

function get_data_end_mean(filename, groupname, probabilities, observer::String)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    h5open(path_to_file, "r") do file
        g_mipt = file[groupname]
        return get_data_end_mean(g_mipt, probabilities, observer)
    end
end
function get_data_end_mean(group::HDF5.Group, probabilities, observer::String)
    out = []
    for p in probabilities
        e_mean = read_attribute(group["p = $p"][observer], "end_mean")
        push!(out, e_mean)
    end
    return out
end

function get_attributes_string(filename, groupname; kwargs...)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    h5open(path_to_file, "r") do file
        g_mipt = file[groupname]
        return get_attributes_string(g_mipt; kwargs...)
    end
end
function get_attributes_string(group::HDF5.Group; addnewline::Bool = true, newlinecharlimit::Integer = 50)
    out = ""
    newlinesadded = 0
    for attr_key in keys(attributes(group))
        value = string(read_attribute(group, attr_key))
        stringToAdd = "$attr_key = $value,"
        if addnewline && length(out * stringToAdd) > newlinecharlimit * (newlinesadded + 1)
            out *= "\n"
            newlinesadded += 1;
        end
        out *= stringToAdd * " "
    end
    return out[1:end-2]
end

function get_attribute(filename, groupname, attribute)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    h5open(path_to_file, "r") do file
        g_mipt = file[groupname]
        return read_attribute(g_mipt, attribute)
    end
end