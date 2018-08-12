from utils.cleanup_alignments import remove_low_mapq_reads
from utils.cleanup_alignments import remove_low_cellcount_reads

# ------------------------- #
# Remove low quality reads

rule remove_low_mapq_reads:
    """Remove reads with low mapping quality"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.bam")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.bam")
    params: minmapq = lambda wc: wc.minmapq
    wildcard_constraints:
       minmapq='\d+'
    run: 
       remove_low_mapq_reads(input[0], output[0], int(params.minmapq))

# ------------------------- #
# Sort the split reads

rule sort_split_reads:
    """Sort split reads"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.bam")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.sorted.bam")
    wildcard_constraints:
       minmapq='\d+'
    shell: "samtools sort {input} -o {output}"

# ------------------------- #
# Index the reads

rule index_minmapq_reads:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.sorted.bam")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.sorted.bam.bai")
    wildcard_constraints:
        minmapq='\d+'
    shell: "samtools index {input}"


# ------------------------- #
# Index the reads

rule index_deduplicate_reads:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.bam")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.bam.bai")
    wildcard_constraints:
        minmapq='\d+'
    shell: "samtools index {input}"

# ------------------------- #
# Deduplicate the split reads by barcode

rule deduplicate_split_reads_by_barcode:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.sorted.bam"), \
           join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.sorted.bam.bai")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.bam")
    wildcard_constraints:
       minmapq='\d+'
    run:
        deduplicate_reads(input[0], output[0])

INPUT_ALL.append(expand(rules.deduplicate_split_reads_by_barcode.output, 
                        reference=config['reference'], 
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq']))

# ------------------------- #
# Remove barcodes with low read counts

rule remove_low_fragmentcount_barcodes:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.bam")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam")
    params: mincounts = lambda wc: wc.mincounts
    wildcard_constraints:
        mincounts='\d+', minmapq='\d+'
    run:
        remove_low_cellcount_reads(input[0], output[0], int(params.mincounts))

INPUT_ALL.append(expand(rules.remove_low_fragmentcount_barcodes.output, 
                        reference=config['reference'], 
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode']))

# ------------------------- #
# Index the reads

rule index_deduplicate_countfiltered_reads:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam.bai")
    wildcard_constraints:
        mincounts='\d+', minmapq='\d+'
    shell: "samtools index {input}"


rule create_bigwig:
    """Create bigwig of alignment"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bw")
    wildcard_constraints:
        mincounts='\d+', minmapq='\d+'
    threads: 10
    shell:
        "bamCoverage -b {input} -o {output} -p {threads}"

INPUT_ALL.append(expand(rules.create_bigwig.output, 
                        reference=config['reference'], 
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode']))
