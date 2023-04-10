# TemplateGen

TemplateGen is a snakebids-based [BidsApp](https://bids-apps.neuroimaging.io/) that flexibly combines images across modalities to create custom template spaces. It can output both participant-specific transforms into the custom space, and cohort wide template outputs in a format compatible with [TemplateFlow](https://www.templateflow.org/). 

This workflow uses greedy instead of ANTS, for the sake of efficiency, in fact, the registrations as compared to `ants_build_template` seem to be more accurate (though that of course is likely just due to parameter selection, not an inherent limitation of ANTS).. There is ~20-30x speedup for a single pairwise registration compared to ANTS, making template-generation for hundreds of subjects a job that can be completed in around a day with modest resources (32cores). The 4-core 16gb memory greedy registration jobs take <30mins.

## Installation

TemplateGen has a few installation options suitable for different environments:

### Singularity

...

### Docker

...

### pip

Pip installations are the most flexible, and offer the best opportunity for parallelization on cluster environments, but are also the most involved, as 3rd party requirements will not be automatically installed. See below for different strategies on dealing with this.

First, be sure that python 3.7 or greater is installed on your system:

```bash
python --version
```

Then start by creating a virtualenv:

```bash
python -m venv .venv
source .venv/bin/activate
```

And install via pip:

```bash
pip install templategen
```

Alternatively, you can get a safe, user-wide installation using [pipx](https://pypa.github.io/pipx/) (pipx installation instructions found [here](https://pypa.github.io/pipx/installation)):

```bash
pipx install templategen
```

#### 3rd party requirement options for pip installation

##### Singularity

With singularity installed on your system, TemplateGen can take care of 3rd party installations by itself. The singularity executable must be available on your `$PATH` (try `singularity --version` on the command line), and your computer must have an active internet connection so that the singularity containers can be downloaded.

This option works well for clusters and compute environments, which typically have singularity already installed.

##### Manual Installation

The required software can also be manually installed. TemplateGen depends on [ANTs v2.3.4](https://github.com/ANTsX/ANTs/releases/tag/v2.3.4) and [itk-SNAP v4.0](http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.SNAP4). In particular, you must have the following commands available on your `$PATH`:

* From ANTs:
    * `AverageImages`
    * `MultiplyImages`
    * `AverageAffineTransformNoRigid`
    * `antsApplyTransforms`
    * `ResampleImageBySpacing`
* From itk-SNAP
    * `greedy`
    * `c3d_-_affine_tool`


## Usage



## Running with a cluster profile (e.g. cc-slurm)
This option will be the most # of jobs (`num_iters * num_subjects * num_cohorts`), but will maximize parallelization on a cluster
```
snakemake  --profile cc-slurm 
```

## Running with cluster profile, but grouping registration jobs together in chunks of 8
This will reduce the number of jobs by factor of 8 by grouping the 8 4core jobs to fill the 32-core nodes

HOWEVER -- if you are using the --group-components option to group the registrations in chunks of 8, you *MUST* run each iteration separately, using the `--config run_iter=#` option, where # is the iteration to run -- this seems to be a bug/limitation of the group-components option when a job has recursively-defined rules. (open issue as of Oct 3 2020: https://github.com/snakemake/snakemake/issues/656)
```
snakemake --config run_iter=2 --profile cc-slurm --group-components reg=8 composite=100
```
Note: composite=100 is for the composing warps to another ref space (e.g. mni) to combine 100 subjects in a single job (each is relatively quick and disk-bound)..


## Running as a single job
This is the most frugal for resources (32-cores max), but is simple in that only a single job is spawned.. This is the method to use if you are running this on a single machine. Four iterations of building a single 1mm template with 100 subjects takes < 24hrs on a 32core system.

regularSubmit -j Fat snakemake --use-singularity  -j32

## Running each cohort in parallel with single jobs and the --nolock flag
This runs each cohort separately with single-node jobs.. This is a happy medium for running on a cluster while not needing hundreds (or thousands) of short jobs.. 

WARNING: make sure your cohorts are mutually-exclusive when using this method, as it is running snakemake in parallel on the same directory with the --nolock option -- if cohorts are not mutually-exclusive, you can still use this method, but only after all the pre-processing (e.g. T2/T1 registration, is completed)
```
for cohort in young middle1 middle2 old1 old2; do  regularSubmit -j Fat snakemake --use-singularity -j32 --nolock --config run_cohort=$cohort; done
```

