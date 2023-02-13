using TOML

function package_version_number()
    path_to_project = joinpath(@__DIR__, "..", "..", "Project.toml")
    return TOML.parsefile(path_to_project)["version"]
end