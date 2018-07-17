from utils.split_reads import augment_alignment_by_barcode_from_name

rule augment_by_barcodes_from_name:
    """Split reads by barcodes"""
    input: join(OUT_DIR, '{reference}', '{sample}.bam')
    output: temp(join(OUT_DIR, "{reference}", "{sample}.barcoded.bam"))
    run:
      augment_alignment_by_barcode_from_name(input[0], output[0])
