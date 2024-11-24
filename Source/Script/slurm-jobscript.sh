snakemake --use-conda --jobs 100 \
          --cluster-config Config/cluster.yaml \
          --cluster "sbatch --partition={cluster.partition} --time={cluster.time} --cpus-per-task={cluster.cpus} --mem={cluster.mem} --output={cluster.output} --error={cluster.error}"
