

#### loads GTEX RSEM isoform level results, filtered for just HK1

function loadgtexgct(metafile="gtex_meta.tsv.gz", tpmfile="gtex_hk1.tsv.gz"; projdir=getprojectdir())

    meta = CSV.read(joinpath(projdir, "data", metafile), DataFrame)

    io = open(joinpath(projdir, "data", tpmfile)) |> GzipDecompressorStream
    
    verline = readline(io)
    
    rows, cols = parse.(Int, split(readline(io)))
    header = split(readline(io))
    
    M = zeros(rows, cols)
    
    genes = String[]
    trans = String[]
    
    for i = 1:rows
        line = readline(io)
        fields = split(line)
        push!(trans, fields[1])
        push!(genes, fields[2])
        
        for j = 1:cols
           M[i, j] = parse(Float64, fields[j+2])
        end
        
    end
    
    close(io)
    
    df = [DataFrame([trans, genes], header[1:2]) DataFrame(M, header[3:end])]
    sdf = stack(df, names(df)[3:end], variable_name=:SampleID, value_name=:TPM)

    
    tpm = innerjoin(meta, sdf, on=:SampleID)  ### isoform level
    genetpm = combine(groupby(tpm, [:SampleID, :Tissue, :gene_id]), :TPM => sum => :TPM) ## gene level

    tissues = sort(unique(genetpm.Tissue))
    
    tissueorder = Dict(t => i for (i, t) in enumerate(tissues));
    mediantissues = sort(combine(groupby(genetpm, :Tissue), :TPM => median => :Median), :Median, rev=true).Tissue
    mediantissueorder = Dict(t => i for (i, t) in enumerate(mediantissues))
    genetpm.TissueOrder = getindex.(Ref(tissueorder), genetpm.Tissue);
    genetpm.MedianOrder = getindex.(Ref(mediantissueorder), genetpm.Tissue);

    tpm.TissueOrder = getindex.(Ref(tissueorder), tpm.Tissue);
    tpm.MedianOrder = getindex.(Ref(mediantissueorder), tpm.Tissue);



    (tpm=tpm, genetpm=genetpm)
end

function gtexisoformheatmap(tpm; c=:blues, num_iso=3, kwargs...)

    tmm = @chain tpm begin
        groupby([:Tissue, :transcript_id])
        combine(:TPM => median => :TPM)
        unstack(:Tissue, :transcript_id, :TPM)
        coalesce.(_, 0)
    end

    M = Matrix(tmm[!, 2:end])
    sir = sortperm(maximum(M, dims=2) |> vec, rev=true)
    sic = sortperm(maximum(M, dims=1) |> vec, rev=true)
    
    ni = 1:num_iso
    heatmap(identity.(M[sir, sic[ni]]), yflip=true, c=c, yticks=(1:size(M, 1), tmm.Tissue[sir]), xticks=(ni, names(tmm)[2:end][sic[ni]]), xrotation=45, 
        size=(500, 800), grid=false, left_margin=40mm, fontfamily="Helvetica", ylabel="Tissue", xlabel="HK1 Isoform", title="Top $(num_iso) most expressed HK1 Isoforms" ; kwargs...)

end

function gtexboxplot(genetpm; rev=false, orderfield=:MedianOrder)
    plot(grid=false)
    tissues = sort(unique(gtex.genetpm[!, [:Tissue, orderfield]]), orderfield, rev=rev).Tissue

    cc = cgrad(:viridis, length(tissues), categorical=true)
    order = genetpm[!, orderfield]
    order = ifelse(rev, maximum(order) .- order .+ 1, order)

    boxplot!(order, genetpm.TPM, group=genetpm[!, orderfield], c=hcat(cc...), marker=(stroke(0), 0.5, 2), size=(1000,600),  lab="")
    plot!(xticks=(1:length(tissues), tissues), xrotation=-45, bottom_margin=50mm, right_margin=10mm, fontfamily="helvetica", ylabel="TPM", title="GTEx HK1 Expression", left_margin=5mm, xlabel="Tissue")
end

function plot_hk1_isoforms(tpm, genecoords; num_iso=3, fs=10)
    tmm = @chain tpm begin
        groupby([:Tissue, :transcript_id])
        combine(:TPM => median => :TPM)
        unstack(:Tissue, :transcript_id, :TPM)
        coalesce.(_, 0)
    end

    M = Matrix(tmm[!, 2:end])
    sic = sortperm(maximum(M, dims=1) |> vec, rev=true)
    topisos = names(tmm)[sic[1:num_iso] .+ 1]
    
    topiso_ind = Int[]
    for t in topisos
        tn = replace(t, r"\.[0-9]$" => "")       
        fi = findall(occursin.(tn, genecoords.TranscriptID))
        @assert length(fi) == 1
        push!(topiso_ind, first(fi))
    end

    p = plot()

    si = sortperm(topiso_ind, by = x -> genecoords.stop[x] - genecoords.start[x], rev=true)


    yos = 0:(num_iso-1)
    
    cc = cgrad(:viridis, num_iso+1, categorical=true)
    cc = [:steelblue, :orange, :black]
    for k in 1:num_iso
        yo = yos[k]
        i = topiso_ind[si[k]]
        plotgene!(genecoords.TranscriptID[i], genecoords, y=1 - yo , off=1, label=genecoords.TranscriptID[i], fs=fs, pos=:left, c=cc[k])
    end

    xl = Int.(round.(xlims()))
    xlt = string("chr10:", first(xl), "-", last(xl))
    plot!(size=(975, 150), bottom_margin=5mm, yticks=false, yaxis=false, left_margin=18mm, ylims=(-3, 3), xlabel=xlt, title="Top Expressed HK1 isoforms", fontfamily="helvetica")
    
end