snakemake -s main.snake --jobs 4 --local-cores 10 --resources mem_mb=1000 \
 --cluster "qsub -V -l h_vmem={resources.mem_mb}m -l h_stack=300m -pe smp {threads} -N {rule} " --reason --rerun-incomplete
