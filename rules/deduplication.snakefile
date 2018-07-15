
# ------------------------- #
# Sort the split reads

rule sort_split_reads:
    """Sort split reads"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.bam")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.sorted.bam")
    shell: "samtools sort {input} -o {output}"

# ------------------------- #
# Deduplicate the split reads by barcode

rule deduplicate_split_reads_by_barcode:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.sorted.bam")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam")
    run:
        deduplicate_reads(input[0], output[0])

INPUT_ALL.append(expand(rules.deduplicate_split_reads_by_barcode.output, reference=config['reference'], sample=config['samples'].keys()))

# ------------------------- #
# Index the reads

rule index_deduplicate_split_reads:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam.bai")
    shell: "samtools index {input}"