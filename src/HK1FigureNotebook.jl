using DataFrames, DataFramesMeta, CSV, CodecZlib
using Plots, StatsPlots, Measures
using BioSequences, FASTX, GenomicFeatures, GenomicIntersections, GenomeFragments
using MotifScanner
using ProgressMeter
using Statistics

include("mutations.jl")
include("plots.jl")
include("motifs.jl")
include("regions.jl")
include("genomics.jl")
include("conformation.jl")
include("expression.jl")

function showwide(table)
    c = ENV["COLUMNS"]
    ENV["COLUMNS"] = "10000"
    display(table)
    ENV["COLUMNS"] = c;
    nothing;
end

function showwl(lines=200, nc=10000)
    
    table -> begin
        l = ENV["LINES"]
        c = ENV["COLUMNS"]
        ENV["LINES"] = string(lines)
        ENV["COLUMNS"] = string(nc)
        display(table)
        ENV["LINES"] = l;
        ENV["COLUMNS"] = c;
        nothing;
    end
end

function getprojectdir()
    d = pwd()
    if basename(d) == "notebooks"
        return dirname(d)
    else
        return d
    end
end