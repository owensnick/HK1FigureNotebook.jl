
using DataFrames, DataFramesMeta, CSV, CodecZlib
using Plots, StatsPlots, Measures
using BioSequences, FASTX, GenomicFeatures, GenomicIntersections, GenomeFragments
using MotifScanner
using ProgressMeter
using Statistics, StatsBase


## Some packages are loaded within included files below
include("mutations.jl")
include("plots.jl")
include("motifs.jl")
include("regions.jl")
include("genomics.jl")
include("conformation.jl")
include("expression.jl")


"""
    loaddata(loadalignments=true, loadhic=true; coolerfile="")

    Loads data necessary for analysis and plotting. This repo contains all source data necessary to regenerate the figures calculated from original alignment files and cooler file.
    Should you wish to recalculate source files, then include alignment files in data/chip_atac_datasets.tsv (.bin format see https://github.com/owensnick/GenomeFragments.jl) and
    the cooler file supplied as an argument to this function. Then set loadalignments=true and loadhic=true.

    
    


"""
function loaddata(loadalignments=false, loadhic=false; coolerfile="")
    
    ### mutations and region
    mutations = loadmutations()
    region = loadregregion() 
    annotateregion!(mutations, region)
    distinctmutations, refseq, altseqs = distinctmuts(mutations, region);
    hk1coords = hk1_locus_coords()
    gtf = loadgtf()
    genecoords = transcriptcoords(gtf);

    ### chip/atac datasets
    datasets = loaddatasets();
    
    
    datagroup = combine(groupby(datasets, [:Study, :DataType, :Cell, :Target]), nrow => :count, :Index => Ref => :Inds);
    datagroup.Label = @with datagroup string.(:Study, "_", :DataType, "_", :Cell, "_", :Target);
    sort!(datagroup, [:Study, :Target, :Cell]);
    
    
    
    if loadalignments
        FM = load_frag_matrix.(datasets.SignalFile);
    else
        FM = nothing
    end
    
    ### hicdata
    hicdata = loadhicdata();
    
    
    ### Expression datasets
    gtex = loadgtexgct();
    betadiffstage = load_beta_diff_geusz_xie();
    betadiffsc = load_betadiff_sc();
    islet_scrnaseq = loadisletscrna();
    betadifflate = load_betadiff_scrnaseq();
    
    
    ### cooler
    if loadhic

        if !isfile(replace(coolerfile, r"::resolutions\/[0-9]*$" => "")) ## strip off cooler resolution
            @warn "Cooler file $coolerfile not found. Proceeding without loading"
            isletcooler = nothing
        else
            isletcooler = loadcoolerfile(coolerfile)
        end

    else
        isletcooler = nothing
    end
    hicdata = quantisletloops(hicdata, isletcooler) ## quantify HIC data in islet loops
    
    mutations, region, distinctmutations, refseq, altseqs, hk1coords, genecoords, datasets, datagroup, FM, hicdata, gtex, betadiffstage, betadiffsc, islet_scrnaseq, betadifflate, isletcooler
end


#### Displaying wide tables, primary use in jupyter notebooks
function showwide(table)
    c = ENV["COLUMNS"]
    ENV["COLUMNS"] = "10000"
    display(table)
    ENV["COLUMNS"] = c;
    nothing;
end


#### Displaying wide tables and long tables, returns a function that will do displaying
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


### function to get base dir of project whether working in src or notebooks
function getprojectdir()
    d = pwd()
    if (basename(d) == "notebooks") || (basename(d) == "src")
        return dirname(d)
    else 
        return d ## assume working in project dir
    end
end