
### useful plotting functions

### savedisplay save svg and display current figure

function savedisplay(label, fmt=".svg", projdir=getprojectdir())
    file = joinpath(projdir, "figures", string(label, fmt))
    mkpath(dirname(file))
    savefig(file)
    p = plot!()
    display(p)
end

######## mutation plots

function plotmutations(region, mutations)
    mutfam = mutations[!, Not([:Patient, :Inheritance])] |> unique
    mutgroups = combine(groupby(mutfam[mutfam.Contained, :], [:RelReg]), nrow => :Count, :type => Ref => :type, :Alt => Ref => :Alt)
    
    p = plotseq(region.metadata)

    acgt = MotifScanner.motif_letter_data(weight="normal")

    for row in eachrow(mutgroups)
        for i = 1:row.Count
            if row.type[i] == "var"
                plotletter!([row.RelReg[1] .- 0.5, -i*0.275], [1.0, 0.25], acgt[first(row.Alt[i])], c=:red) ### hard coded sizes to match MotifScanner, TODO: remove these
            end

            if row.type[i] == "del"
                plot!(row.RelReg, -0.25*ones(length(row.RelReg)), c=:red, lab="")
                plot!(row.RelReg[1]*[1, 1], [-0.1, -0.25], c=:red, lab="")
                plot!(row.RelReg[end]*[1, 1], [-0.1, -0.25], c=:red, lab="")
                annotate!([(mean(row.RelReg), -0.5, text("del", font("helvetica oblique", 9, :red)))])
            end
        end
    end
    p
end


"""
    idconvert(x)

    Function to trim lead 711086.. off snp Ids
"""
function idconvert(x)
    fields = split(x, "_")
    string(replace(fields[1], "711086" => ""), ":", fields[2], ">", fields[3])
end

function disruptmotifscatter(motsel ; kwargs...)

    theme(:wong2)
    
    @with motsel scatter(:RefPrMax, :AltPrMax, group=:MotifFamily, marker=(stroke(0)))#, marker=^(:auto))
    plot!([0, 1], [0, 1], c=:black, lab="", leg=:outertopright, xlabel="Motif Score Ref.", ylabel="Motif Score Mut")
    xlims!(0.5, 0.8)
    ylims!(-.125, 0.8)
    yticks!(-0.0:.2:0.8)
    annotate!([(motsel.RefPrMax[i] .+ 0.01, motsel.AltPrMax[i], text(motsel.ID[i] |> idconvert, "helvetica", :left, 6)) for i = 1:2:size(motsel, 1)])
    annotate!([(motsel.RefPrMax[i] .- 0.01, motsel.AltPrMax[i], text(motsel.ID[i] |> idconvert, "helvetica", :right, 6)) for i = 2:2:size(motsel, 1)])
    ps = plot!(size=(275, 225), grid=false, title="Top Ranked Motifs\ndisrupted by mutation", fontfamily="helvetica", titlefont=font(10, "helvetica"); kwargs...)
    theme(:wong)
    ps

end

function disruptedmotif_seqlogo(region, mutations, motsel)
    p = plotmutations(region, mutations)
    plot!(size=(600, 200), ylims=(-2, 3.5), grid=false, yaxis=false, yticks=false)
    
    for row in eachrow(motsel)
        motif = motifs[findfirst(motifmeta.MotifName .== row.MotifName)]
        xo = row.RefSeqStart-1
        yo = 0.5 #+ ifelse(row.MotifName .== "Hic1", 2, 0)
        seqlogo!(motif, rc=row.RefStrand == "-", xo=xo, yo=yo)
        annotate!((xo + size(motif.pwm, 2)/2, yo + 2.5, text(row.MotifFamily, font("helvetica", :top, :center, 8))))
    end
    xl = (71108635, 71108685+12)
    plot!(xticks=(1:25:101, string.( (1:25:101) .+ leftposition(region) .- 1)), xlims = xl .- leftposition(region) .+ 1, xlabel=string("chr10:", first(xl), "-", last(xl)))
    plot!(bottom_margin=5mm, fontfamily="helvetica")
    p
end



##### heatmap averaging

function avh(H::Vector{T}, δ) where {T}
    n = length(H)

    w = cld(n, δ)
    AH = zeros(Float64, w)
    te = zeros(Int, w)
    for i = 1:n
        wi = cld(i, δ)
        te[wi] += 1
        AH[wi] += H[i]
    end
    if te[end] < 0.75δ
        AH[end] += AH[end-1]
        te[end] += te[end-1]
    end
    AH./te
end


############### plotting gene models

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
function plotgene!(trans, exoncoords;  y=1, off=1, h = 0.5, label="", fs=10, pos=:left, shiftstrand=true, c=:black, kwargs...)
    
    g = @subset(exoncoords, :TranscriptID .== trans)
    p = plot!(size=(1200, 100) ; kwargs...)
    
    if shiftstrand    
        if g.strand[1] == "-"
            # y = y - 1
            # c = :red
        end
    end

    
    for (s, e) in g.coords[1]
       plot!(rectangle(e - s + 1, h, s - off + 1, y), c=c, linecolor=c, lab="")
    end
    plot!([g.start[1], g.stop[1]] .- off .+ 1, y .+ h/2 .+ [0, 0], c=c, lab="")
    if !isempty(label)
        if pos == :left
            annotate!([(g.start[1] - 2000, y + h/2, text(label, font(:right, fs, "helvetica")))])
        else
            ### corrections to yoff
            yc = 0
            if g.GeneName[1] == "HK1"
                yc = -1
            elseif g.GeneName[1] .== "TYSND1"
                yc = 1.5
            end
            annotate!([((g.start[1] + g.stop[1])/2, y-0.5 + yc, text(label,font(:center, fs, "helvetica")))])
        end
    end
    p
end

overlaps(startA, stopA, interval) = !isempty(intersect(startA:stopA, interval))

locs(s, e) = s:e
function plotgenemodels(loc, exoncoords; off=1, fs=10)

    maxmodels = @chain exoncoords begin
        @subset(overlaps.(:start, :stop, Ref(loc)))
        groupby(:GeneName)
        combine(df -> df[argmax(df.stop - df.start), :])
    end
    
    maxmodels.ovg = overlappinglocations(fill("", size(maxmodels, 1)), locs.(maxmodels.start .- 1000, maxmodels.stop .+ 1000))
    maxmodels = transform(groupby(maxmodels, :ovg), :ovg => (x -> 1:length(x)) => :yoff)
    maxmodels[maxmodels.GeneName .== "TACR2", :yoff] .+= 1
    p  = plot()
    
    for (g, t, yo) in zip(maxmodels.GeneName, maxmodels.TranscriptID, maxmodels.yoff)
        plotgene!(t, exoncoords, y=1 - yo + 1, off=off, label=g, fs=fs, pos=:bottom)
    end
    if "ENST00000359426.7_2" ∈ exoncoords.TranscriptID
        plotgene!("ENST00000359426.7_2", exoncoords, y=0, off=off, label="", fs=fs, pos=:bottom) ## extra isoform for hk1 shortest and highested expressed in islets
        #plotgene!("ENST00000298649.8_2", exoncoords, y=-1, off=off, label="", fs=fs, pos=:bottom) ##  second shortest, there is TF binding right on the promoter of this, but does not appear to be accessible or expressed
    end
               
    xlims!(first(loc), last(loc))
    p
end


