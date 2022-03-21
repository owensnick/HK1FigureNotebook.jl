### functions to define useful regions of interest and load gene models

### differeing regions of interest
### coords are 1-based inclusive
function hk1_locus_coords()
    chrom = "chr10"
    mutations = 71108610:71108717

    genelocus = (71029740 - 4000):(71161640+4000)
    
    roi = 71107650:71110242
    criticalregion = 71108645:71108686
    criticalregion_4k = (71108645 - 2000):(71108686 + 2000)
    roi_200k = (first(roi) - 200_000):(last(roi) + 200_000)
    Dict("mutations"         => (chrom, mutations),
         "roi"               => (chrom, roi),
         "roi_200k"          => (chrom, roi_200k),
         "genelocus"         => (chrom, genelocus),
         "criticalregion"    => (chrom, criticalregion),
         "criticalregion_4k" => (chrom, criticalregion_4k),
         "widecontactview"   => (chrom, 70257650:71960242),
         "ins_contactview"   => ("chr11", 1330974:3032599),
         "sst_contactview"   => ("chr3", 186536699:188238191),
         "gcg_contactview"   => ("chr2", 162149392:163858914))
end



### load gtf file
function loadgtf(;file="gencode.v36lift37.annotation.gois.gtf.gz", projdir=getprojectdir(), fp = joinpath(projdir, "data", file))

    gtf = CSV.read(fp, DataFrame, comment="#", header=[:chrom, :source, :feature, :start, :stop, :dot, :strand, :frame, :description]);
    
    gtf = @subset(gtf, :feature .== "exon")
    dd = descriptiondict.(gtf.description)
    gtf[!, :GeneID] = getindex.(dd, "gene_id")
    gtf[!, :GeneName] = getindex.(dd, "gene_name")
    gtf[!, :Gene] = string.(gtf.GeneID, "|", gtf.GeneName);
    gtf[!, :TranscriptID] = getindex.(dd, "transcript_id")
    gtf[!, :TranscriptType] = getindex.(dd, "transcript_type")

    @subset(gtf, :TranscriptType .== "protein_coding")
end


### turn gtf file description field into a dictionary
function descriptiondict(des)
    fields = split(des, ";", keepempty=false)
    rc = [split(f, keepempty=false) for f in fields]
    Dict(id => replace(rec, "\"" => "") for (id, rec) in rc)
end


transcriptcoords(gtf) = combine(groupby(gtf, [:GeneName, :TranscriptID, :TranscriptType, :strand]), :start => minimum => :start, :stop => maximum => :stop, nrow => :count, [:start, :stop] => ((a, b) -> [collect(zip(a, b))]) => :coords)
