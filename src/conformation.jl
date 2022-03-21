## HiC data + loading cooler code

include("cooler.jl")


function loadhicdata(loadfuns = [load_islet_hic]) 

    dfs = [f() for f in loadfuns]
    ns = names(first(dfs))
    for df in dfs[2:end]
        ns = intersect(ns, names(df))
    end
    mapreduce(df -> df[!, ns], vcat, dfs)
end

function load_islet_hic(file="islet_loops_TSTFF938730.bedpe.gz", projdir=getprojectdir())
    filepath = joinpath(projdir, "data", file)
    isletloops = CSV.read(filepath, DataFrame)
    rename!(isletloops, "#chrA" => "chrA")
    isletloops[!, :Study] .= "IsletHiC"
    isletloops
end
