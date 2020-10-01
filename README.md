# Snakemake workflow: `ants_build_template_camcan`

Large workflow, so needs to be run in steps..

Processing between ANTs runs can be done on a single 32core node, using regularInteractive or regularSubmit.

The ANTS runs should use a slurm job for each run, (`--profile cc-slurm`) to maximize parallelization.. if resources are low or time is not a big issue, can also run this on single 32core node, running in sequence.. 


e.g: running the first iteration:

Run pre-processing (everything up to the first ants registration):
```
regularSubmit -j ShortFat snakemake -j32 --use-singularity all_iter1 --omit-from reg_to_template
```

Run ants registration for iteration 1 (all cohorts, submitted in parallel):
```
snakemake --profile cc-slurm all_iter1  --until reg_to_template
```


