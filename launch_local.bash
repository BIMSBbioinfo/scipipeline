export OMP_NUM_THREADS=1
LOGDIR=/scratch/AG_Akalin/wkopp/alison_sciatac_output/logs
snakemake -s main.snake --jobs 10 --local-cores 10 \
 --cluster "qsub -v OMP_NUM_THREADS -V -l h_vmem={resources.mem_mb}m -l h_stack=10m -pe smp {threads} -N {rule}.\${{JOB_ID}} -o ${LOGDIR}/rule_{rule}.o\${{JOB_ID}} -e  ${LOGDIR}/rule_{rule}.e\${{JOB_ID}} " --reason --rerun-incomplete
