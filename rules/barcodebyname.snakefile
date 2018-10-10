from utils.split_reads import augment_alignment_by_barcode_from_name

rule augment_by_barcodes_from_name:
    """Split reads by barcodes"""
    input: 
      barbam = join(OUT_DIR, '{reference}', '{sample}.cleanchrom.bam'),
      reftable = config['barcodes']['sheet']
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.bam")
    resources:
      mem_mb=500
    run:
      augment_alignment_by_barcode_from_name(input.barbam, 
                                             output[0],
                                             input.reftable)
