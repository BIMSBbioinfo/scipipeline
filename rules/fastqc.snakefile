
rule quality_control:
    """Quality control with fastqc"""
    input: get_mapping_inputs
    output: join(OUT_DIR, 'fastqc', '{sample}')
    shell:
      "mkdir -p {output}; fastqc {input} -o {output}"

INPUT_ALL.append(expand(rules.quality_control.output, sample=config['samples'].keys()))
