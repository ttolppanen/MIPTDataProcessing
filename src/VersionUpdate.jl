# using HDF5

export update_0_1_to_0_2

function update_0_1_to_0_2(filename)
    path_to_file = joinpath(getdatapath(), filename * ".h5")
    h5open(path_to_file, "cw") do file
        file_version = split(read_attribute(file, "version"), ".")
        if file_version[1] == "0" && file_version[2] == "1"
            for g in file
                attributes(g)["save_before_effect"] = false
            end
            delete_attribute(file, "version")
            write_attribute(file, "version", "0.2.0")
        end
    end
end