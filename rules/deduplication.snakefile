from utils.cleanup_alignments import remove_low_mapq_reads
from utils.cleanup_alignments import remove_low_cellcount_reads
from utils.split_reads import deduplicate_reads

# ------------------------- #
# Remove low quality reads

rule remove_low_mapq_reads:
    """Remove reads with low mapping quality"""
    input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.bam")
    output: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.bam"),
            join(OUT_DIR, "{sample}", "{reference}", 'report', "summary_removed_minmapq{minmapq}.tsv")
    params: minmapq = lambda wc: wc.minmapq
    wildcard_constraints:
       minmapq='\d+'
    resources:
        mem_mb=1000
    run:
       #"samtools view -bq {params.mainmapq}"
       remove_low_mapq_reads(input[0], output[0], int(params.minmapq), output[1])

# ------------------------- #
# Sort the split reads

rule sort_split_reads:
    """Sort split reads"""
    input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.bam")
    output: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.sorted.bam")
    wildcard_constraints:
       minmapq='\d+'
    resources:
        mem_mb=5000
    shell: "samtools sort {input} -o {output}"

# ------------------------- #
# Index the reads

rule index_minmapq_reads:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.sorted.bam")
    output: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.sorted.bam.bai")
    wildcard_constraints:
        minmapq='\d+'
    resources:
        mem_mb=2000
    shell: "samtools index {input}"


# ------------------------- #
# Index the reads

rule index_deduplicate_reads:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.bam")
    output: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.bam.bai")
    wildcard_constraints:
        minmapq='\d+'
    resources:
        mem_mb=2000
    shell: "samtools index {input}"

# ------------------------- #
# Deduplicate the split reads by barcode

rule deduplicate_split_reads_by_barcode:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.sorted.bam"), \
           join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.sorted.bam.bai")
    output: outbam=join(OUT_DIR, '{sample}', "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.bam"), \
       summary=join(OUT_DIR, '{sample}', "{reference}", "report", "summary_custom_markdup.minmap{minmapq}.txt")
    wildcard_constraints:
       minmapq='\d+'
    resources:
        mem_mb=5000
    run:
        deduplicate_reads(input[0], output.outbam, output.summary, True)
#rule deduplicate_split_reads_by_barcode:
#    """Deduplicate split reads"""
#    input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.sorted.bam"), \
#           join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.sorted.bam.bai")
#    output: outbam=join(OUT_DIR, '{sample}', "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.bam"), \
#       summary=join(OUT_DIR, '{sample}', "{reference}", "report", "summary_picard_markdup.minmap{minmapq}.txt")
#    params: picard=config['picard_jarpath']
#    threads: 10
#    log: join(LOG_DIR, 'picard_markdup_{sample}_{reference}_minmap{minmapq}.log')
#    wildcard_constraints:
#       minmapq='\d+'
#    resources:
#        mem_mb=5000
#    shell:
#        "java -Xms1000m -Xmx{resources.mem_mb}m -XX:ParallelGCThreads={threads} -jar {params.picard} MarkDuplicates I={input[0]} O={output.outbam} M={output.summary} BARCODE_TAG=RG REMOVE_DUPLICATES=true 2> {log}"

INPUT_ALL.append(expand(rules.deduplicate_split_reads_by_barcode.output,
                        reference=config['reference'],
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq']))


if config['compute_library_complexity']:
    # ------------------------- #
    # Measure library complexity before deduplication
    
    rule library_complexity_before_dedup:
        """Deduplicate split reads"""
        input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.sorted.bam")
        output: join(OUT_DIR, "{sample}", "{reference}", "report", "summary_library_complexity_beforededup.minmap{minmapq}.txt")
        params: picard=config['picard_jarpath']
        threads: 10
        resources:
            mem_mb=5000
        log: join(LOG_DIR, 'picard_estlibcompl_beforededup_{sample}_{reference}_minmapq{minmapq}.log')
        wildcard_constraints:
           minmapq='\d+'
        shell:
            "java -Xms100m -Xmx400m -XX:ParallelGCThreads={threads} -jar {params.picard} EstimateLibraryComplexity I={input} O={output} 2> {log}"
    
    INPUT_ALL.append(expand(rules.library_complexity_before_dedup.output,
                            reference=config['reference'],
                            sample=samples.Name.tolist(),
                            minmapq=config['min_mapq']))
    
    # ------------------------- #
    # Measure library complexity after deduplication
    
    rule library_complexity_after_dedup:
        """Deduplicate split reads"""
        input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.bam")
        output: join(OUT_DIR, "{sample}", "{reference}", "report", "summary_library_complexity_afterdedup.minmap{minmapq}.txt")
        params: picard=config['picard_jarpath']
        log: join(LOG_DIR, 'picard_estlibcompl_afterdedup_{sample}_{reference}_minmapq{minmapq}.log')
        threads: 10
        wildcard_constraints:
           minmapq='\d+'
        resources:
            mem_mb=5000
        shell:
            "java -Xms100m -Xmx400m -XX:ParallelGCThreads={threads} -jar {params.picard} EstimateLibraryComplexity I={input} O={output} 2> {log}"
    
    INPUT_ALL.append(expand(rules.library_complexity_after_dedup.output,
                            reference=config['reference'],
                            sample=samples.Name.tolist(),
                            minmapq=config['min_mapq']))

# ------------------------- #
# Remove barcodes with low read counts

rule remove_low_fragmentcount_barcodes:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.bam")
    output: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam"), \
            join(OUT_DIR, "{sample}", "{reference}", 'report', "summary_mincount_filter.minmapq{minmapq}.mincount{mincounts}.tsv")
    params: mincounts = lambda wc: wc.mincounts
    resources:
        mem_mb=2000
    wildcard_constraints:
        mincounts='\d+', minmapq='\d+'
    run:
        remove_low_cellcount_reads(input[0], output[0], int(params.mincounts), output[1])

INPUT_ALL.append(expand(rules.remove_low_fragmentcount_barcodes.output,
                        reference=config['reference'],
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode']))

# ------------------------- #
# Index the reads

rule index_deduplicate_countfiltered_reads:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam")
    output: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam.bai")
    wildcard_constraints:
        mincounts='\d+', minmapq='\d+'
    resources:
        mem_mb=2000
    shell: "samtools index {input}"


rule create_bigwig:
    """Create bigwig of alignment"""
    input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam"),
           join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam.bai")
    output: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bw")
    wildcard_constraints:
        mincounts='\d+', minmapq='\d+'
    threads: 5
    resources:
        mem_mb=10000
    shell:
        "bamCoverage -b {input[0]} -o {output} -p {threads}"

INPUT_ALL.append(expand(rules.create_bigwig.output,
                        reference=config['reference'],
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode']))
