
rule quality_control_trimmed:
    """Quality control with fastqc"""
    input: get_mapping_inputs
    output: directory(join(OUT_DIR, '{sample}', 'fastqc_trimmed'))
    resources:
      mem_mb=10000
    threads: 1
    shell:
      "mkdir -p {output}; fastqc {input} -o {output}"

INPUT_ALL.append(expand(rules.quality_control_trimmed.output, sample=samples.Name.tolist()))

rule quality_control_raw:
    """Quality control with fastqc"""
    input: get_trim_inputs
    output: directory(join(OUT_DIR, '{sample}', 'fastqc_raw'))
    resources:
      mem_mb=10000
    threads: 1
    shell:
      "mkdir -p {output}; fastqc {input} -o {output}"

INPUT_ALL.append(expand(rules.quality_control_raw.output, sample=samples.Name.tolist()))
