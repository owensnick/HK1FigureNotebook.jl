### functions for plotting gene expression data

### Load GP libraries for psuedotime
using GaussianProcesses
using Distributions

### include 
include("gtex.jl")


#### islet sc rnaseq
function loadisletscrna(metafile="scrnaseq_human_islets_meta.tsv", countfile="scrnaseq_human_islets_normcounts.tsv.gz"; projdir=getprojectdir())
    meta = CSV.read(joinpath(projdir, "data", metafile), DataFrame)
    normcounts = CSV.read(joinpath(projdir, "data", countfile), DataFrame)
    (; meta, normcounts)
end

### plot boxplot of islet sc data
function scboxplot(gene, data ; kwargs...)    
    ind = findfirst(occursin.(gene, data.normcounts.Gene)) 
    plot(title=string("GSE101207\n", gene), leg=false, xlabel="Islet Cell Type", ylabel="Norm. counts", fontfamily="helvetica" , xrotation=45; kwargs...)
    y = [v for v in data.normcounts[ind, data.meta.Label]]

    celllabels = Dict("pp" => "Pancreatic polypeptide", "delta" => "Delta", "psc" => "Pancreatic stellate", "duct" => "Duct")
    celltype = get.(Ref(celllabels), data.meta.CellType, data.meta.CellType)
    boxplot!(celltype, y, outliers=false, c=:grey, linecolor=:black)
    dotplot!(celltype, y, c=:steelblue, marker=(stroke(:white)))
end

function scboxplotgenes(genes, data ; kwargs...)
    phs = [scboxplot(g, islet_scrnaseq, size=(1200, 400)) for g in genes]
    for p in phs[2:end]
        plot!(p, ylabel="", left_margin=-2mm)
    end
    plot!(phs[1], ylabel="Norm. Counts", left_margin=5mm)
    
    p = plot(phs..., layout=(1, length(phs)), size=(1100, 300),  bottom_margin=10mm, top_margin=5mm, fontfamily="helvetica", titlefont=font(12, "helvetica"); kwargs...)
end



### beta cell maturation scRNA-seq
function load_betadiff_scrnaseq(metafile="scrnaseq_balboa_xin_cell_stage_meta.tsv", countfile="scrnaseq_balboa_xin_cell_stage_exp.tsv.gz"; projdir=getprojectdir())
    meta = CSV.read(joinpath(projdir, "data", metafile), DataFrame)
    normcounts = CSV.read(joinpath(projdir, "data", countfile), DataFrame)

    (; meta, normcounts)
end

## bar plot
function barbetadiff(gene, data ; kwargs...)
    ind = findfirst(data.normcounts.GeneName .== gene)
    meta = @subset(data.meta, :cell_type .∈ Ref(["Adult_Alpha", "Adult_Beta"]))
    y = [data.normcounts[ind, l] for l in meta.Label]

    stages = ["S5",  "S6",  "S7",  "M1",  "M3",  "M6",  "A"]
    stageorder = Dict(s => i for (i, s) in enumerate(stages))
    cc = cgrad(:cividis, 9, categorical=true)[2:end]

    groupedbar(meta.cell_type, y, group=getindex.(Ref(stageorder), meta.developmental_stage), lab=permutedims(stages), c=permutedims(cc), title=string("scRNA-seq endocrine cell maturation\n$(gene)"),
        ylabel="Norm. Counts", xlabel="Cell Type", fontfamily="helvetica"; kwargs...)
end
    
######### beta cell differentiation Pseudotime

function load_betadiff_sc(; projdir=getprojectdir())
    
    labelsfile = joinpath(projdir, "data", "singlecell_betadiff_labels.csv")
    tpmfile    = joinpath(projdir, "data", "singlecell_betadiff_tpm.csv.gz")

    data = CSV.read(tpmfile, DataFrame)
    rename!(data, ["Gene" ; string.("t", names(data)[2:end])])

    labels = CSV.read(labelsfile, DataFrame)
    labels.mid = (labels.start .+ labels.stop)/2;
    labels
    

    (data=data, labels=labels)

end

### beta cell differnentiation by stages, not currently used in figures
function load_beta_diff_geusz_xie(; projdir=getprojectdir())

    metafile = joinpath(projdir, "data", "rnaseq_betadiff_geusz_xie_meta.tsv")
    tpmfile  = joinpath(projdir, "data", "rnaseq_betadiff_geusz_xie_tpm.tsv.gz")

    meta = CSV.read(metafile, DataFrame)
    tpm = CSV.read(tpmfile, DataFrame)

    meta.Cell = replace.(meta.Cell, Ref("endocrine" => "EN"))

    (meta=meta, tpm=tpm)
end

function plotgenestage(gene, data ; kwargs...)

    ind = findall(occursin.(gene, data.tpm.Gene))
    
    if isempty(ind)
        error("Gene $gene not found")
    else
        plotgenestage(ind, data; kwargs...) 
    end
end
function plotgenestage(ind::Vector{T}, data; maxnorm=false, vert=false, cc=1:length(ind), kwargs...) where {T}
    
    meta = data.meta
    tpm = data.tpm
    tdf = unique(sort(meta[!, [:CellOrder, :Cell]], :CellOrder))
    p = plot()
    
    
    if  length(ind) > 1
        ofv = range(.1, -.1, length=length(ind))
        s = 2
    else
        ofv = [0]
        s = 1
    end

    for (k, i) in enumerate(ind)
        x = [tpm[i, l] for l in meta.Label]
        if maxnorm
            x ./= maximum(x)
        end
        mg = combine(groupby(DataFrame(CellOrder=meta.CellOrder, X=x), :CellOrder), :X => mean => :Mean, :X => std => :Std)        
        
        if vert

            scatter!(x, meta.CellOrder .+ meta.Offset/s .+ ofv[k], yticks=(tdf.CellOrder, tdf.Cell), c=cc[k], linecolor=cc[k], lab=tpm.GeneName[i], leg=:outertopright, marker=(1, stroke(0));  kwargs...)
            plot!(mg.Mean, mg.CellOrder .+ ofv[k], marker=:hline, xerror=mg.Std, c=cc[k], lab="")
        else
            scatter!(meta.CellOrder .+ meta.Offset, x, xticks=(tdf.CellOrder, tdf.Cell), c=cc[k], lab="", leg=:outertopright, marker=(stroke(0));  kwargs...)
            plot!(mg.CellOrder, mg.Mean, marker=:hline, yerror=mg.Std, lab=tpm.GeneName[i], c=cc[k])
        end

    end
    p
end


########## GaussianProcess Regression for sc beta diff psuedotime

datatrans(x, α=100, β = 1) = log(α*x + β)
invtrans(y, α=100, β=1) = (exp(y)-β)/α
sqd(x) = ifelse(x > 0, x*x, 0);
function gp_reg(t, y, st; initparams = (1.0, 0.0, 4.0), datatrans=datatrans, invtrans=invtrans, f_prior = Normal(1.4, 4.0), ℓ_prior = Normal(2., 2.0), n_prior=Normal(1.0, .75), label="")
    σf, ℓ, σn = initparams

    kern = Mat52Iso(ℓ, σf)
    set_priors!(kern, [ℓ_prior, f_prior])
    gp = GP(t, datatrans.(y), MeanZero(), kern, σn)
    set_priors!(gp.logNoise, [n_prior] )
    optimize!(gp)
    μ, v = predict_y(gp, st)
    sf  = invtrans.(μ)
    cil = invtrans.(μ .- 1.96*sqrt.(v))
    ciu = invtrans.(μ .+ 1.96*sqrt.(v))

    (gp=gp, μ=μ, v=v, sf=sf, cil=cil, ciu=ciu, t=t, y=y, st=st, label=label)
end

function gpgene(gene, data, dx=0.5)
    x = collect(@subset(data, :Gene .== "PseudoTime")[1, 2:end]) 
    y = collect(@subset(data, :Gene .== gene)[1, 2:end]); 
    y[end] = mean(y[end-3:end-1]) ## final datapoint appears outlying, replace with mean of previous 3
    st = minimum(x):dx:maximum(x)
    gp_reg(x, y, st, label=gene)
end

function plotgpr(gpr ; c=:steelblue)
    plot()
    scatter!(gpr.t, gpr.y, c=c, lab="", marker=(stroke(0), 2))
    plot!(gpr.st, gpr.sf, ribbon=[gpr.sf - gpr.cil gpr.ciu - gpr.sf], lab=gpr.label, c=c)
end

function plotgpr!(gpr ; c=:steelblue, maxnorm=false, vert=true)

    if maxnorm
        m = maximum(gpr.sf)
    else
        m = 1
    end

    if vert
        
        scatter!(gpr.y/m, gpr.t, c=c, lab="", marker=(stroke(0), 2))
        plot!(gpr.sf/m, gpr.st, lab="", c=c)

        plot!(Shape([gpr.cil/m ; reverse(gpr.ciu/m)], [gpr.st ; reverse(gpr.st)]), fillalpha=0.2 ,line=stroke(0), c=c, lab=gpr.label)
    else
        

        scatter!(gpr.t, gpr.y/m, c=c, lab="", marker=(stroke(0), 2))
        plot!(gpr.st, gpr.sf/m, ribbon=[gpr.sf/m - gpr.cil/m gpr.ciu/m - gpr.sf/m], lab=gpr.label, c=c)
    end
end


function plotgprgenes(genes, betadiffsc ; kwargs...)
    phs = Plots.Plot[]
    for (i, g) in enumerate(genes)
        p = plot()
        plotgpr!(gpgene(g, betadiffsc.data), c=i, vert=false)
        p = xticks!(betadiffsc.labels.mid, betadiffsc.labels.cell)
        plot!(title=g)
        push!(phs, p)
    end
    plot!(phs[1], ylabel="Norm counts", left_margin=5mm)
    for i in 2:length(phs)
        plot!(phs[i], yformatter=y->"", left_margin=-3mm) 
    end
    p = plot(phs..., layout=(1, length(phs)), size=((200*length(phs)) + 70, 195), link=:y, leg=false, grid=false, right_margin=-2mm, framestyle=:box, fontfamily="helvetica"; kwargs...)

end


function gp_hk1_gck(betadiffsc)
    plot(size=(325, 175), xticks=0:2, ylabel="Expression (au)", leg=:top, grid=false, xlabel="Differentiation Pseudotime")
    plotgpr!(gpgene("HK1", betadiffsc.data), maxnorm=true, vert=false)
    plotgpr!(gpgene("GCK", betadiffsc.data), c=:orange, maxnorm=true, vert=false)
    xticks!(betadiffsc.labels.mid, betadiffsc.labels.cell)
    plot!(fontfamily="helvetica", yticks=0:2)
end