
### Functions for loading/anntoting mutations/variants

using FASTX

### load mutations

function loadmutations(; file="hk1mutations.tsv", projdir=getprojectdir())
    filepath = joinpath(projdir, "data", file)
    df = CSV.read(filepath, DataFrame)

    df.Ref = [ismissing(s) ? dna"" : LongDNASeq(s) for s in df.Ref]
    df.Alt = [ismissing(s) ? dna"" : LongDNASeq(s) for s in df.Alt]

    df.PatientGroup = Int.(floor.(df.Patient))
    df[!, :ID] = string.(df.start, "_", df.Ref, "_", df.Alt)

    df
end

function parsechromloc(s)
    fields = split(replace(s, "," => ""), r"[:-]")
    fields[1], parse(Int, fields[2]):parse(Int, fields[3])
end


function loadregregion(;file="hk1_reg_region.fa",  projdir=getprojectdir())
    filepath = joinpath(projdir, "data", file)
    reader = open(FASTA.Reader, filepath)
    record = first(collect(reader))
    seq = sequence(record)
    chromloc = replace(first(split(description(record))), "range=" => "")
    chrom, loc = parsechromloc(chromloc)
    seq
    Interval(chrom, loc, '+', seq)
end


function annotateregion!(mutations, region)

    regc = leftposition(region):rightposition(region)
    mutreg = [s:e for (s, e) in zip(mutations.start, mutations.stop)]
    relreg = [mr .- leftposition(region) .+ 1 for mr in mutreg]
    mutations[!, :RelReg] = relreg
    mutations[!, :Contained] = length.(intersect.(Ref(regc), mutreg)) .== length.(mutreg)

    ### assert that ref alleles match sequence
    refseqs = [region.metadata[r] for r in mutations[mutations.Contained, :RelReg]]
    @assert refseqs == mutations.Ref[mutations.Contained]
    mutations

end

buildmutation(refreg, kind, ref, alt, contained) = (coords=refreg, kind=kind, ref=ref, alt=alt, contained=contained)

function buildaltseq(seq, mutation)
    if !mutation.contained
        return missing, missing, missing, missing
    end
    leftind = 1:(first(mutation.coords) - 1)
    refrightind = (last(mutation.coords) + 1):length(seq)
    if mutation.kind == "var"
        altind  = mutation.coords
        altrightind = refrightind
    elseif mutation.kind == "del"
        altind = first(mutation.coords) .+ (0:-1)
        altrightind = first(mutation.coords):(length(seq) - length(mutation.coords))
    else
        error("Mutation type: $mutation.kind not recognised")
    end

    
    altseq = seq[leftind]*mutation.alt*seq[refrightind]
    leftind, altind, refrightind, altrightind, altseq
end

function distinctmuts(mutations, region)
    distinctmutations = combine(groupby(@subset(mutations, :Contained), names(mutations, Not([:Patient, :Mutation, :Inheritance]))), nrow=>:count)
    
    muts = buildmutation.(distinctmutations.RelReg, distinctmutations.type, distinctmutations.Ref, distinctmutations.Alt, distinctmutations.Contained);
        
    refseq = region.metadata
    altseqindex = buildaltseq.(Ref(refseq), muts)
    distinctmutations.LeftInd = first.(altseqindex)   
    distinctmutations.AltInd = getindex.(altseqindex, 2)
    distinctmutations.RefRightInd = getindex.(altseqindex, 3)
    distinctmutations.AltRightInd = getindex.(altseqindex, 4)
    
    altseqs = last.(altseqindex)
        
    distinctmutations, refseq, altseqs
end