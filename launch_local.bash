export OMP_NUM_THREADS=1
snakemake -s main.snake --jobs 10 --local-cores 10 \
 --cluster "qsub -v OMP_NUM_THREADS -V -l h_vmem={resources.mem_mb}m -l h_stack=10m -pe smp {threads} -N {rule}.${JOB_ID} " --reason --rerun-incomplete
