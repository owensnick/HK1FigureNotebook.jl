## Functions for scanning for motifs that may be disrupted by variants

### build scanning sequence for a given motif and mutation described by indexes to the left and to the right exclusive of mutation
### this allows scanning of all positions in which motif intersects with at least one bp of mutated sequence
function motifscanseq(seq, leftind, rightind, motif)
    n = size(motif.pwm, 2)
    start = max(last(leftind) - n + 1 + 1, 1)
    stop  = min(first(rightind) + n - 1 - 1, length(seq))
    seq[start:stop], start
end


### find maximal forward or reverse match of a sequence
function scanmax(seq, motif)
    fs, rs = scanmotif(seq, motif.pbg)
    fm, fi = findmax(fs)
    rm, ri = findmax(rs)
    n = size(motif.pbg, 2)
    if fm > rm
        maxscore = fm
        start = fi
        stop = fi + n - 1
        strand = "+"
    
    else
        maxscore = rm
        start = ri
        stop  = ri + n - 1
        strand = "-"
    end
    maxscore, start, stop, strand
end

### scan all mutations with all motifs
function motifscannall(refseq, altseqs, distinctmutations, motifs)

    dfs = DataFrame[]
    for (as, row) in zip(altseqs, eachrow(distinctmutations))
        df = scanmots(refseq, as, row.LeftInd, row.RefRightInd, row.AltRightInd, motifs)
        df[!, :ID] .= row.ID
        push!(dfs, df)
    end
    reduce(vcat, dfs)
end


### scan a given mutation with all motifs
function scanmots(refseq, altseq, leftind, refrightind, altrightind, motifs)
    df = DataFrame(MotifName=String[], MotifID=String[], MotSeqStart=Int[], RefMaxScore=Float64[], RefStart=Int[], RefStop=Int[], RefStrand=String[], AltMaxScore=Float64[], AltStart=Int[], AltStop=Int[], AltStrand=String[],
        RefPrMax=Float64[], AltPrMax=Float64[], LR_RefAlt=Float64[], PR_RefAlt=Float64[])
    
    for m in motifs
        maxscore = sum(maximum(m.pbg, dims=1))
        refmseq, refstart = motifscanseq(refseq, leftind, refrightind, m)
        altmseq, altstart = motifscanseq(altseq, leftind, altrightind, m)
        refres = scanmax(refmseq, m)
        altres = scanmax(altmseq, m)

        ref_prmax = first(refres)/maxscore
        alt_prmax = first(altres)/maxscore
        lr_refalt = first(refres) - first(altres)
        pr_refalt = ref_prmax - alt_prmax

        
        push!(df, (m.name, m.id, refstart, refres..., altres..., ref_prmax, alt_prmax, lr_refalt, pr_refalt))
    end
    df.RefSeqStart = df.MotSeqStart + df.RefStart .- 1
    df.RefSeqStop  = df.MotSeqStart + df.RefStop .- 1

    df
end


### load motifs from JASPAR transfac file
function loadtransfac(file="JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac.txt.gz", w = 1e-2, background=MotifScanner.background_human_zero_markov(), projdir=getprojectdir())

    filepath = joinpath(projdir, "data", file)
    io = open(filepath) |> GzipDecompressorStream
    mode = 1
    
    motifname = ""
    motifid = ""
    currentmotif = Vector{Vector{Float64}}()
    pwms = Vector{Matrix{Float64}}()
    nsites = Float64[]
    
    motifdata = DataFrame(MotifName=String[], MotifID=String[], Property=String[], Value=String[])
    
    for line in eachline(io)
        
        fields = split(line)
        
        if mode == 1
            if fields[1] == "AC"
                mode = 2
                currentmotif = Vector{Vector{Float64}}()
            end
        elseif mode == 2
            if fields[1] == "DE"
                motifid = fields[2]
                motifname = fields[3]
                if startswith(motifname, motifid) ## JASPAR 2022 release appears to have altered the format of JASPAR 2020, remove the leading motifid
                    motifname = replace(motifname, string(motifid, ".") => "")
                end
            elseif fields[1] == "PO"
                @assert fields[2] == "A"
                @assert fields[3] == "C"
                @assert fields[4] == "G"
                @assert fields[5] == "T"
                mode = 3
            elseif fields[1] == "CC"
                ind = findfirst(':', line)
                
                cfields = [line[4:(ind-1)], line[(ind+1):end]]
                
                if length(cfields) != 2
                    @show line
                    @show cfields
                    error("")
                end
                push!(motifdata, (motifname, motifid, cfields...))
            elseif fields[1] == "//"
                mode = 1
            end
        elseif mode == 3
            if fields[1] == "XX"
                mode = 2
                pwm = reduce(hcat, currentmotif)
                v = sum(pwm, dims=1)
                push!(nsites, mean(v))
                pwm ./= sum(pwm, dims=1)
                pwm .+= w
                pwm ./= sum(pwm, dims=1)
                push!(pwms, pwm)
            else
                v = parse.(Float64, fields[2:end])
                push!(currentmotif, v)
            end
        end
    end
    close(io)
    motifmeta = unstack(motifdata, :Property, :Value)

    
    motifmeta[!, :NumSites] = nsites
    motifmeta[!, :Length] = size.(pwms, 2)

    @assert iszero(mapreduce(c -> sum(ismissing, c), +, eachcol(motifmeta)))
    dropmissing!(motifmeta)
    
    pbgs = [log2.(pwm) .- log2.(background) for pwm in pwms]

    motifdata = [(name=motifname, id=id, pwm=pwm, pbg=pbg) for (motifname, id, pwm, pbg) in zip(motifmeta.MotifName, motifmeta.MotifID, pwms, pbgs)]
    motifmeta.MotifFam = motiffamily.(motifmeta.MotifName);
    motifmeta, motifdata

end


### motif family from motif tf name, crude but does the job
function motiffamily(motifname)
    mf = uppercase.(motifname)
    mf = replace(mf, r"\([A-Za-z0-9.]*\)$" => "")
    
    mf = replace(mf, r"[0-9]*$" => "")
    mf = replace(mf, r"-$" => "")
    mf = replace(mf, r"FOX[A-Z]" => "FOX")
end

function moteq(mA, mB)
    (mA.name != mB.name) && return false
    (mA.id != mB.id) && return false
    
    (maximum(mA.pwm - mB.pwm) > 1e-5) && return false
    return true
end

overlap(grp, starts, stops) = overlappinglocations(fill(grp, length(starts)), [s:e for (s, e) in zip(starts, stops)])
function filtermatches(motscan; minscore=0.6, maxscore=Inf)

    motmatch = @chain motscan begin
        @subset(minscore .<= :RefPrMax .< maxscore)                                   # filter by norm motif score
        groupby(:ID)                                                                  # group by mutation id
        transform(:PR_RefAlt => (x -> sortperm(sortperm(x, rev=true))) => :Rank)      # Rank all motif matches (grouped by ids)
        sort([:MotifFamily, :ID, :Rank])                                              # sort
        transform([:MotifFamily, :RefSeqStart, :RefSeqStop] => overlap => :Group)     # Assign grouping label for overlapping motifs
        sort([:ID, :Rank])                                                            # sort 
    end
    motmatch
end

function motif_rank1(motmatch)
    motsel = @chain motmatch begin
        @subset(:Rank .== 1)                                                          # rank of motif for each mutation
        groupby(:Group)                                                               # multiple rank 1 motifs overlap, group by those that overlap
        combine(df -> df[argmax(df.PR_RefAlt), :])                                    # select the motif with largest disruption score
        sort(:ID)
    end
    motsel
end

function motifmatchtable(motmatch_tier1, motmatch_tier2, region)

    motmatch_tier1[!, :Tier] .= 1
    motmatch_tier2[!, :Tier] .= 2

    motmatchall = [ motmatch_tier1 ; motmatch_tier2]
    
    motmatchall.MutSeqStart = motmatchall.MotSeqStart .+ motmatchall.AltStart .- 1
    motmatchall.MutSeqStop = motmatchall.MotSeqStart .+ motmatchall.AltStop .- 1
    motmatchall.RefSeqMotifStart = leftposition(region) .+ motmatchall.RefSeqStart .- 1
    motmatchall.RefSeqMotifStop  = leftposition(region) .+ motmatchall.RefSeqStop .- 1
    motmatchall.MutSeqMotifStart = leftposition(region) .+ motmatchall.MutSeqStart .- 1
    motmatchall.MutSeqMotifStop  = leftposition(region) .+ motmatchall.MutSeqStop .- 1
    motmatchall[!, :chrom] .= "chr10"
    motmatchall[!, :Mutation] = [string.("g.", f[1], ifelse(f[3] == "", string("del", f[2]), string(f[2], ">", f[3]))) for f in split.(motmatchall.ID, "_")] 
    rename!(motmatchall, :RefStrand => :RefSeqStrand, :RefMaxScore => :RefSeqMaxScore, :AltStrand => :MutSeqStrand, :AltMaxScore => :MutSeqMaxScore,
                :RefPrMax => :RefSeqNormMaxScore, :AltPrMax => :AltSeqNormMaxScore, :LR_RefAlt => :LR_RefMut, :PR_RefAlt => :Norm_RefMut)
    mma = motmatchall[!, [:Mutation, :Tier, :MotifFamily, :MotifName, :MotifID, :chrom,
                        :RefSeqMotifStart, :RefSeqMotifStop, :RefSeqStrand, :RefSeqMaxScore, 
                        :MutSeqMotifStart, :MutSeqMotifStop, :MutSeqStrand, :MutSeqMaxScore, 
                        :RefSeqNormMaxScore, :AltSeqNormMaxScore, :LR_RefMut, :Norm_RefMut, :Rank, :Group]] 
    
    mma
end