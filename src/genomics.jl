### large number of functions for loading, analysing and plotting genomic data e.g. ChIP-seq, ATAC-seq etc.

function loaddatasets(file="chip_atac_datasets.tsv", projdir=getprojectdir())
    meta = CSV.read(joinpath(projdir, "data", file), DataFrame)

    sort!(meta, [:Study, :Target, :Cell, :Rep])

    meta = @subset(meta, :Target .!= "Input");
    meta.Index = 1:size(meta, 1)
    meta
end

#### Functions to generate pileups of Chip/ATAC at loci and cache results
pile(chrom, loc, FM::FragMatrixSingle{T}, fl=120) where{T}  = 1e+6*fragregion(chrom, loc, "+", FM, inc_meta_frag!, fl, 4, one)/totalfrags(FM)
pile(chrom, loc, FM::FragMatrixPair{T}, minfrag=0, maxfrag=1000) where{T}  = 1e+6*fragregion(chrom, loc, FM, inc_meta_frag!, minfrag, maxfrag, 4, one)/totalfrags(FM)

function loadpiles(file, labels)
    df = CSV.read(file, DataFrame)
    @assert names(df) == labels
    collect(eachcol(df))
end

function piles_saveload(chrom, loc, inds, labels, FM::Nothing, file; saveload=true)
    !isfile(file) && error("Error: precalcuated $file not found, load alignments to regenerate")
    if saveload == false
        @warn "Loading precalculated signal tracks set to false, but alignments not loaded, loading precalculated signals"
    end
    loadpiles(file, labels)
end
function piles_saveload(chrom, loc, inds, labels, FM, file; saveload=true)
    saveload && isfile(file) && return loadpiles(file, labels)
    P = @showprogress map(ind -> mapreduce(i -> vec(pile(chrom, loc, FM[i])), +, ind), inds)

    if saveload
        mkpath(dirname(file))
        CSV.write(file, DataFrame(P, labels), delim='\t', compress=true)
    end
    P
end
function piles(chrom, loc, inds, grouplabel, labels, FM; saveload=true, projdir=getprojectdir())
    file = joinpath(projdir, "results", "piles", string("pile_", grouplabel, "_", chrom, "_", first(loc), "_", last(loc), ".tsv.gz"))
    piles_saveload(chrom, loc, inds, labels, FM, file, saveload=saveload)
end
piles(regionlabel, grouplabel, datagroup::DataFrame, regioncoords, FM ; saveload=true) = piles(regioncoords[regionlabel]..., datagroup.Inds, grouplabel, datagroup.Label, FM, saveload=saveload)



### generate a vector plots of pile ups from different data sources
function prp(chrom, loc, datagroup, P, regioncoords; dx=50, annot=true, regspace=0, plotrec=true, plotline=false, cc=[], ylm=-1)

    labels = string.(datagroup.Cell, " ", datagroup.Target)
    
    pp = first(loc):dx:last(loc)
    phs = Plots.Plot[]
    for (i, p) in enumerate(P)
        c = isempty(cc) ? :auto : cc[i]

        q = plot()
        if ylm != -1
            yl = (0, ylm)
            plot!(ylims=(yl))
        else
            yl = ylims()
        end

        ### Mark pileup with ROI containing variants
        if plotrec && haskey(regioncoords, "roi")
            roi = last(regioncoords["roi"])
            plot!(rectangle(last(roi) - first(roi), yl[2] - yl[1], first(roi), 0), fillcolor=:lightgrey, linecolor=nothing, lab="", fillalpha=0.5)
        end
        if plotline && haskey(regioncoords, "criticalregion")
            cr = last(regioncoords["criticalregion"])
            vline!([first(cr), last(cr)], c=:black, lab="", ls=:dot)
        end
        plot!(pp, avh(p, dx), lab="", fill=0, c=c)
        plot!(ylims=(yl))
        
       
        annot && annotate!((first(loc), (yl[1] + yl[2])/2, text(labels[i], font(:left, 10, :middle, "helvetica"))))
        if i < length(P)
            plot!(q, xticks=false, bottom_margin=-3Measures.mm)
        end        
        push!(phs, q)
    end
    phs
end


##################### functions to generate different genomics panels plots

####### histones #######
function calc_histonepanel(regionlabel, datagroup, regioncoords, FM, genecoords, saveload=true)
    datahistone = @subset(datagroup, occursin.(r"H3K27ac|H3K27me3|H3K4me1|H3K4me3", :Target), :Study .!= "PPCebola")

    cellorder = Dict("ES" => 1, "DE" => 2, "GT" => 3, "PP1" => 4, "PP2" => 5, "EndoCβH1" => 6, "islets" => 7, "Alpha" => 8, "Beta" => 9, "Exocrine" => 10)
    datahistone.CellOrder = getindex.(Ref(cellorder), datahistone.Cell)
    sort!(datahistone, [:Target, :CellOrder])
    datahistone.GroupIndex = 1:size(datahistone, 1)

    PP = piles(regionlabel, "histone", datahistone, regioncoords, FM, saveload=saveload)
    p = histonepanel(regioncoords[regionlabel]..., datahistone, PP, genecoords, regioncoords, cc=cgrad(:viridis, length(unique(datahistone.Cell)) + 2, categorical=true)[2:end-1])
end

function histonepanel(chrom, loc, datahistone, PH, genecoords, regioncoords; dxs=[500, 500, 1000, 500], cc=cgrad(:viridis, length(unique(datahistone.Cell)), categorical=true))
    histoneorder = ["H3K27ac", "H3K4me1", "H3K27me3"]
    ylms = [1.0, 1.5, 1.0, 1]
    phs = Vector{Plots.Plot}[]
    for (h, ylm, dx) in zip(histoneorder, ylms, dxs)
        dh = @subset(datahistone, :Target .== h)
        ps = prp(chrom, loc, dh, PH[dh.GroupIndex], regioncoords, dx=dx, annot=false, regspace=2000, cc=cc)
        
        for (p, l) in zip(ps, dh.Cell)
            plot!(p, ylims=[0, ylm], yticks=0:1)
            if l == "ES"
                plot!(title=h)
            end
            if h == first(histoneorder)
                plot!(p, ylabel=l, left_margin=5Measures.mm)
            end
        end

        push!(phs, ps)
    end

    mp = maximum(length.(phs))


    ### fill in grid with empty plots
    for ps in phs
        for i = (length(ps)+1):mp
            plot!(ps[end], xticks=false, bottom_margin=-3Measures.mm)
            p = plot(xlims=xlims(ps[end]), yaxis=false, yticks=false, framestyle=:none)
            push!(ps, p)
        end
    end


    for ps in phs
        plot!(ps[end], xticks=false, bottom_margin=-3Measures.mm)
        p = plotgenemodels(loc, @subset(genecoords, :GeneName .== "HK1"))
        plot!(ylims=(-1, 2), yticks=false, yaxis=false, xticks=false, xaxis=false, xlabel=string(chrom, ":", first(loc), "-", last(loc)))
        push!(ps, p)
    end

    
    plot(permutedims(reduce(hcat, phs))..., layout=(length(phs[1]), length(histoneorder)), link=:x, fontfamily="helvetica", xlims=(first(loc), last(loc)), size=(1000, 400), grid=false)

end

####### pancreatic progenitors #######

function calc_progentior_panel(regionlabel, datagroup, regioncoords, FM, genecoords)
    dataprog = @subset(datagroup, .!occursin.(r"^H[23]", :Target), :DataType .== "ChIP", :Study .!= "EndoCChIP", :Cell .∉ Ref(["DE", "ES", "GT", "islets"]), :Target .∉ Ref(["CTCF"]), 
    (:Study .!= "PPCebola") .| (:Target .== "TEAD1") .| (:Cell .== "LiverBud"))

    dataprog = [sort(@subset(dataprog, .!occursin.(r"FOX", :Target)), :Target) ; sort(@subset(dataprog, occursin.(r"FOX", :Target)), [:Target, order(:Cell, rev=true)])]
    dataprog.GroupIndex = 1:size(dataprog, 1)

    PP = piles(regionlabel, "pancprog", dataprog, regioncoords, FM)
    pancpropanel(regioncoords[regionlabel]..., dataprog, PP, genecoords, regioncoords)
    p = plot!(size=(492, 400))# xlims=(71065952-5000,71065952+5000))

end

function pancpropanel(chrom, loc, dataprog, PP, genecoords, hk1coords)

    yl = (0, 3)
    theme(:wong)
    cc= getindex.(Ref(Dict(u => i for (i, u) in enumerate(unique(dataprog.Target)))), dataprog.Target)
    cc = replace(cc, 4 => maximum(cc)+1)
    phs = prp(chrom, loc, dataprog, PP, hk1coords, dx=100, annot=false, ylm=yl[2], cc=cc)
    labels = string.(dataprog.Target, " ", ifelse.(occursin.("PP", dataprog.Cell), "PP", "Liver Bud"))
    
    for (i, p) in enumerate(phs)
        annotate!(p, [(first(loc), (yl[1] + yl[2])/2, text(labels[i], font(:left, 10, :middle, "helvetica")))])
        plot!(p, yticks=[0, 3])
    end

    plot!(phs[end], xticks=false, bottom_margin=-3Measures.mm)

    p = plotgenemodels(loc, @subset(genecoords, :GeneName .== "HK1"))
    plot!(ylims=(-1, 2), yticks=false, yaxis=false, xticks=false,  xaxis=false, xlabel=string(chrom, ":", first(loc), "-", last(loc)))
    push!(phs, p)

    plot(phs..., layout=(length(phs), 1), size=(700, 400), grid=false, link=:x, xlims=(first(loc), last(loc)))
end


####### islets #######
function calcisletpanel(regionlabel, datagroup, regioncoords, FM, genecoords)

    dataislets = @subset(datagroup, .!occursin.(r"^H[23]", :Target),
                                    :Study .∉ Ref(["DiffChip", "DiffChipATAC", "EndoCChIP", "PPCebola", "EndoCATAC"]),
                                    (:Study .!= "PancSNATAC") .| (:Cell .== "beta_1"))
    
    PP = piles(regionlabel, "islet_chip", dataislets, regioncoords, FM)

    ## sort islet data by region in ROI descending

    roi_iv = last(regioncoords["roi"]) .- first(last(regioncoords[regionlabel])) .+ 1
    atac_ind = dataislets.DataType .== "ATAC"

    si = sortperm(map(p -> maximum(p[roi_iv]), PP[.!atac_ind]), rev=true)
    PP = [PP[.!atac_ind][si] ; PP[atac_ind]]
    dataislets = [dataislets[.!atac_ind, :][si, :] ; dataislets[atac_ind, :]]

    
    isletpanel(regioncoords[regionlabel]..., dataislets, PP, genecoords, regioncoords)
    plot!(size=(480, 350))
end


function isletpanel(chrom, loc, datagroup, PP, genecoords, hk1coords)

    yl = (0, 10)
    theme(:wong)
    cc= getindex.(Ref(Dict(u => i for (i, u) in enumerate(unique(datagroup.Target)))), datagroup.Target)
    cc = replace(cc, 4 => maximum(cc)+1)
    phs = prp(chrom, loc, datagroup, PP, hk1coords, dx=100, annot=false, ylm=yl[2], cc=cc)
    labels = ifelse.(datagroup.DataType .== "ATAC", "snATAC beta-cells", datagroup.Target)
    
    for (i, p) in enumerate(phs)

        if datagroup.DataType[i] == "ATAC"
            ylt = (0, 3)
        else
            ylt = yl  
        end
        plot!(p, yticks=[ylt...], ylims=ylt)
        annotate!(p, [(first(loc), (ylt[1] + ylt[2])/2, text(labels[i], font(:left, 10, :middle, "helvetica")))])
    end

    plot!(phs[end], xticks=false, bottom_margin=-3Measures.mm)


    p = plotgenemodels(loc, @subset(genecoords, :GeneName .== "HK1"))
    plot!(ylims=(-1, 2), yticks=false, yaxis=false, xticks=false, xaxis=false, xlabel=string(chrom, ":", first(loc), "-", last(loc)))
    push!(phs, p)

    plot(phs..., layout=(length(phs), 1), size=(700, 400), grid=false, link=:x, fontfamily="helvetica")
end

####### snATAC #######
function calc_snatachpanel(regionlabel, datagroup, regioncoords, FM, genecoords)
    datasnatac_hh = @subset(datagroup, :Target .== "ATAC", :Study .== "PancSNATAC", occursin.(r"alpha_1|beta_1|delta_1", :Cell))
    datasnatac_hl = @subset(datagroup, :Target .== "ATAC", :Study .== "PancSNATAC", occursin.(r"alpha_2|beta_2|delta_2", :Cell))

    PPH = piles(regionlabel, "snatac_hormone_high", datasnatac_hh, regioncoords, FM)
    PPL = piles(regionlabel, "snatac_hormone_low",  datasnatac_hl, regioncoords, FM)

    snatacpanel(regioncoords[regionlabel]..., datasnatac_hh, datasnatac_hl, PPH, PPL, genecoords, regioncoords)
end


function snatacpanel(chrom, loc, datasnatac_hh, datasnatac_hl, PPH, PPL, genecoords, hk1coords)

    cc = cgrad(:inferno, 4, categorical=false)
    phs_H = prp(hk1coords["genelocus"]..., datasnatac_hh, PPH, hk1coords, dx=100, ylm=5, cc=cc)
    phs_L = prp(hk1coords["genelocus"]..., datasnatac_hl, PPL, hk1coords, dx=100, ylm=5, cc=cc)

    plot!(phs_H[end], xticks=false, bottom_margin=-3Measures.mm)
    plot!(phs_L[end], xticks=false, bottom_margin=-3Measures.mm)
    p = plotgenemodels(last(hk1coords["genelocus"]), @subset(genecoords, :GeneName .== "HK1"))
    plot!(ylims=(-1, 2), yticks=false, yaxis=false, xticks=false,  xaxis=false)
    
    for (h, l) in zip(phs_L, phs_H)
        plot!(h, yticks=[0, 5]) 
        plot!(l, yticks=[0, 5]) 
    end
    push!(phs_L, p)
    push!(phs_H, deepcopy(p))

    
    plot(mapreduce(collect, vcat, zip(phs_H, phs_L))..., layout=(length(phs_L), 2), size=(1000, 200), grid=false, xlims=(first(loc), last(loc)), link=:x, fontfamily="helvetica")
end

####### atac diff #######
function calcatacdiff_subpanel(regionlabel, datagroup, regioncoords, FM, genecoords)

    data_atac_diff = @subset(datagroup, :DataType .== "ATAC", :Cell .∈ Ref(["ES", "DE", "GT", "PP1", "PP2", "beta_1"]))
    data_atac_diff = sort(data_atac_diff, order(:Cell, by=x -> ifelse(x == "ES", "A", ifelse(x == "EndoC", "Z", x))))

    PP = piles(regionlabel, "atac_diff", data_atac_diff, regioncoords, FM)
    cc = cgrad(:viridis, size(data_atac_diff, 1)+2, categorical=false)

    chrom, loc = regioncoords[regionlabel]
    phs = prp(chrom, loc, data_atac_diff, PP, regioncoords, dx=100, ylm=3.5, cc=cc[2:end-1], annot=false, plotrec=false, plotline=true)
    for (p, l) in zip(phs, data_atac_diff.Cell)
        plot!(p, title=l, yticks=false, left_margin=-1mm)
    end


    plot!(phs[3], xlabel=string(chrom, ":", loc[1], "-", loc[end]))
    plot!(phs[1], yticks=[0, 3], left_margin=5mm)

    p = plot(phs..., layout=(1, length(phs)), size=(325, 100), grid=false, right_margin=-1mm, yaxis=false, bottom_margin=5mm, fontfamily="helvetica",
            tick_direction=:out, xticks=([mean(loc)], ""), titlefont=font(12, "helvetica"))

end

function calc_atacdiff_supp(regionlabel, datagroup, regioncoords, FM, genecoords)

    data_atac_diff = @subset(datagroup, :DataType .== "ATAC", :Cell .∈ Ref(["ES", "DE", "GT", "PP1", "PP2"]))
    data_atac_diff = sort(data_atac_diff, order(:Cell, by=x -> ifelse(x == "ES", "A", ifelse(x == "EndoC", "Z", x))))

    PP = piles(regionlabel, "atac_diff", data_atac_diff, regioncoords, FM)
    cc = cgrad(:viridis, size(data_atac_diff, 1)+1, categorical=false)

    phs = prp(regioncoords[regionlabel]..., data_atac_diff, PP, regioncoords, dx=100, ylm=5, cc=cc)
    plot!(phs[end], xticks=false, bottom_margin=-3Measures.mm)
    for p in phs
        plot!(p, yticks=[0, 5])
    end

    p = plotgenemodels(last(regioncoords[regionlabel]), @subset(genecoords, :GeneName .== "HK1"))
    plot!(ylims=(-1, 2), yticks=false, yaxis=false, xticks=false,  xaxis=false)

    push!(phs, p)
    p = plot(phs..., layout=(length(phs), 1), size=(492, 280), grid=false, link=:x)
end



