export setdatapath
export getdatapath

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