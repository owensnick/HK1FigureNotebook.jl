### Working with cooler files, uses cool https://github.com/open2c/cooler

using PyCall
using DataFrames
using OffsetArrays, CoordinateTransformations, Rotations
using ImageTransformations

loadcooler() = pyimport("cooler")
loadcoolerfile(coolerfile; cooler=loadcooler()) = cooler.Cooler(coolerfile)
loadisletscooler(coolerfile = "C:\\home\\projects\\islets\\hic\\4dn\\islets_DFF064KIG.mcool::resolutions/5000", cooler=loadcooler()) = cooler.Cooler(coolerfile)


### load bin dirs
coolercoords(chrom, loc) = (replace(chrom, "chr" => ""), first(loc), last(loc))
function bindf(chrom, loc, cool; ext=0.5, saveload=true, projdir=getprojectdir())

    dffile =  joinpath(projdir, "results", "cooler", string("isletcooler_bindf_", chrom, "_", first(loc), "_", last(loc), "", ".tsv.gz"))
    if isnothing(cool) || !saveload
        if !isfile(dffile)
            error("Error: cannot load precalculated contact matrix bins $dffile, please load cooler file")
        end
        return CSV.read(dffile, DataFrame)
    end

    dx = Int(round((last(loc) - first(loc))*ext))
    xloc = (first(loc) - dx):(last(loc) + dx)

    bdf = cool.bins().fetch(coolercoords(chrom, xloc))
    cols = Symbol.(bdf.columns)
    data = [[v for v in bdf[c]] for c in cols]
    df = DataFrame(data, cols)
    df.binmid = div.(df.start .+ df.end, 2)
    df

    if saveload
        mkpath(dirname(dffile))
        CSV.write(dffile, df, compress=true)
    end
    df
end


## cooler matrix
function coolmatrix(chrom, loc, cool; balance="KR", ext=0.5, saveload=true, projdir=getprojectdir())

    matfile =  joinpath(projdir, "results", "cooler", string("isletcooler_", chrom, "_", first(loc), "_", last(loc), "", ".tsv.gz"))
    if isnothing(cool) || !saveload
        if !isfile(matfile)
            error("Error: cannot load precalculated contact matrix file $matfile, please load cooler file")
        end
        return Matrix{Float64}(CSV.read(matfile, DataFrame)[!, 2:end])
    end

    dx = Int(round((last(loc) - first(loc))*ext))
    xloc = (first(loc) - dx):(last(loc) + dx)
    CM = cool.matrix(balance=balance).fetch(coolercoords(chrom, xloc))

    if saveload
        mkpath(dirname(matfile))
        CSV.write(matfile, DataFrame(CM, :auto), compress=true)
    end
    CM
end

### rotate matrix to be 2d
function rotcooler(xp, CM ; kwargs...)
    R = warp(CM, recenter(RotMatrix(-π/4), [0, 0]) ; kwargs...)
    n = cld(size(R, 1), 2)
    R = R.parent[(n+1):end, :]
    nxp = range(first(xp), last(xp), length=size(R, 2))
    dx = xp[2] - xp[1]
    nyp = dx*(1:size(R, 1))
    nxp, nyp, R
end

## quantify islet loops
function quantisletloops(hicdata, cool, binsize=5000, ext=2, saveload=true, projdir=getprojectdir())
    datafile = joinpath(projdir, "results", "cooler", string("islet_hicquant.tsv.gz"))
    if isnothing(cool) || !saveload
        if !isfile(datafile)
            error("Error: cannot load precalculated quantification of islet loops $datafile, please load cooler file")
        end
        return CSV.read(datafile, DataFrame)
    end

    dx = binsize*ext
    μ = Float64[]
    @showprogress for row in eachrow(hicdata)
        locA = (row.startA - dx):(row.stopA + dx)
        locB = (row.startB - dx):(row.stopB + dx)
        m = isletcooler.matrix(balance="KR").fetch(coolercoords(row.chrA, locA), coolercoords(row.chrB, locB))
        push!(μ, mean(m))
    end

    hicdata.LoopQuant = μ
    if saveload
        mkpath(dirname(datafile))
        CSV.write(datafile, hicdata, delim='\t', compress=true)
    end
    hicdata
end



# function hk1coolerpanel(chrom, loc, isletcooler, hicdata, regioncoords, ext=0.75)

#     bincoords = bindf(chrom, loc, isletcooler, ext=ext)
#     CM = coolmatrix(chrom, loc, isletcooler, ext=ext)
#     nxp, nyp, R = rotcooler(bincoords.binmid, CM)
#     cl = (0, 2.3139)
#     ph = heatmap(nxp, nyp, log10.(R .+ 1), ylims=(0, 4.0e+5), colorbar=false, c=:inferno, xticks=false, clims=cl) ## clims hardcoded to compare to HK1 scale
#     # cl = zlims()
#     # @show cl
#     yticks!((0:200:400)*1e+3, string.(0:.2:.4))

#     hicloops = @subset(hicdata, :Study .== "IsletHiC", :chrA .== chrom, (:startA .∈ Ref(loc)) .| (:stopB .∈ Ref(loc)) .| (:startB .∈ Ref(loc)) .| (:stopA .∈ Ref(loc))) ## islet loops
#     plottriloops!(hicloops, c=:white, lab="Islet Loops", plotlines=true)

#     plot!(ylabel="Distance (Mb)")

#     vline!([mean(last(regioncoords["roi"]))], c=:white, lab="ROI", ls=:dash, title="Islets HiC")
#     pg = plotgenemodels(loc, genecoords, fs=6)
#     plot!(ylims=(-2, 2), yticks=false, yaxis=false, xticks=true, xlabel=string(chrom, ":", loc[1], "-", loc[2]), bottom_margin=5mm)
#     vline!([mean(last(hk1coords["roi"]))], c=:black, lab="", ls=:dash)
#     pblank = plot([0], [0], zcolor=[NaN], colorbar=true, clims=cl, lab="", framestyle=:none, c=:inferno)
#     lt = @layout [[a{0.65h}  ; b] c{0.05w}]
    
#     p = plot(ph, pg, layout=lt, pblank, size=(900, 290), link=:x, xlims=(first(loc), last(loc)), fmt=:png, background_color_legend=:black, foreground_color_legend=:white, legendfontcolor=:white,
#     left_margin=5mm, fontfamily="helvetica")

# end

### Plot HiC matrix over a region
function coolerpanel(chrom, loc, isletcooler, hicdata, regioncoords;  ext=0.75, plotroi=true)


    bincoords = bindf(chrom, loc, isletcooler, ext=ext)
    CM = coolmatrix(chrom, loc, isletcooler, ext=ext)
    nxp, nyp, R = rotcooler(bincoords.binmid, CM)
    cl=(0, 2.3139)
    ph = heatmap(nxp, nyp, log10.(R .+ 1), ylims=(0, 4.0e+5), colorbar=false, c=:inferno, xticks=false, clims=cl) ## clims hardcoded to compare to HK1 scale
    # cl = zlims()
    # @show cl
    yticks!((0:200:400)*1e+3, string.(0:.2:.4))

    hicloops = @subset(hicdata, :Study .== "IsletHiC", :chrA .== chrom, (:startA .∈ Ref(loc)) .| (:stopB .∈ Ref(loc)) .| (:startB .∈ Ref(loc)) .| (:stopA .∈ Ref(loc))) ## islet loops
    plottriloops!(hicloops, c=:white, lab="Islet Loops", plotlines=true)

    plot!(ylabel="Distance (Mb)")

    plotroi && vline!([mean(last(regioncoords["roi"]))], c=:white, lab="ROI", ls=:dash, title="Islets HiC")
    pg = plotgenemodels(loc, genecoords, fs=6)
    plot!(ylims=(-2, 2), yticks=false, yaxis=false, xticks=true, xlabel=string(chrom, ":", loc[1], "-", loc[2]), bottom_margin=5mm)
    plotroi && vline!([mean(last(hk1coords["roi"]))], c=:black, lab="", ls=:dash)
    pblank = plot([0], [0], zcolor=[NaN], colorbar=true, clims=cl, lab="", framestyle=:none, c=:inferno)
    lt = @layout [[a{0.65h}  ; b] c{0.05w}]
    
    p = plot(ph, pg, layout=lt, pblank, size=(900, 290), link=:x, xlims=(first(loc), last(loc)), fmt=:png, background_color_legend=:black, foreground_color_legend=:white, legendfontcolor=:white,
    left_margin=5mm, fontfamily="helvetica")

end



### plot triangluar loops over HiC Matrix
function plottriloops!(df; w=15000, np=100, c=:white, ss=1, plotlines=true, lab="")
   
    ps = Plots.partialcircle(0, 2π, np, w)
    ofs = ss*w/sqrt(2)
    
    for (k, row) in enumerate(eachrow(df))
        l = (row.startA + row.stopA)/2
        r = (row.startB + row.stopB)/2

        m, s = coordconv(l, r)
        plotlines && plot!([l, m - ofs], [0, s-ofs*sqrt(2)], c=c, lab="")
        plotlines && plot!([m+ofs, r], [s-ofs*sqrt(2), 0], c=c, lab="")
        plot!(first.(ps) .+ m, last.(ps).*sqrt(2) .+ s, c=c, lab=ifelse(k == 1, lab, ""))
    end
    plot!()
end

function coordconv(start, stop) 
    m = (start+stop)/2
    s = (stop - start)/sqrt(2)
    m, s
end