# Snakemake Tutorial with Spack

This is a small repository that shows running the snakemake tutorial using
spack. For the full tutorial, you should see [the snakemake documentation](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html). This tutorial is a subset of that, and we skip over the details
to explain the differences of running natively or with conda and with spack.

## Install Snakemake

For this tutorial, we are currently working off of the [add/spack](https://github.com/snakemake/snakemake/tree/add/spack) 
branch of the snakemake repository. This will eventually be integrated with the
snakemake software, but for now is there. This means that you should currently
clone the branch:

```bash
git clone -b add/spack https://github.com/snakemake/snakemake
cd snakemake
```

Check out a virtual environment:

```bash
python -m venv env
source env/bin/activate
```

And install snakemake locally! This is considered a development install so you
don't muck around with your system or main snakemake.

```bash
$ pip install -e .
```

## Clone the workflow

You can then clone the tutorial here.

```bash
$ git clone https://github.com/snakemake/snakemake-spack-tutorial
cd snakemake-spack-tutorial
```

Take a look around - you'll see an input data directory (data), a script
directory with the Python script to make a plot, and then a folder of spack
environments.

```bash
├── data
├── envs
├── README.md
├── scripts
└── Snakefile
```

## Changes to the Snakefile

The Snakefile is our main configuration file or recipe for defining and
then executing our workflow. If you compare this workfile to the 
[one in the tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#summary),
the only change is that the `conda` directives are now `spack`. Here is an example"

```yaml
rule bwa:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        temp("mapped/{sample}.bam")
    spack:
        "envs/spack-mapping.yml"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

And the environments folder includes yaml files that define spack environments!
Here is the one that we see above:

```yaml
spack:
    specs:
      - bwa
      - samtools
    concretization: together
```

## Running the workflow

Since we need spack, this means that it should be on your path before starting.
You can check as follows:

```bash
$ which spack
```

If you need help installing spack, see [here](https://spack.readthedocs.io/en/latest/).
The spack executable is located in the bin of the repository that you clone.
We then can run the workflow!

```bash
$ snakemake --cores 1 --rerun-incomplete --use-spack
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	all
	2	bwa
	1	call
	2	sort
	1	stats
	7

[Fri May 14 10:56:54 2021]
rule bwa:
    input: data/genome.fa, data/samples/B.fastq
    output: mapped/B.bam
    jobid: 5
    wildcards: sample=B

[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 25000 sequences (2525000 bp)...
[M::mem_process_seqs] Processed 25000 reads in 0.709 CPU sec, 0.709 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 1 data/genome.fa data/samples/B.fastq
[main] Real time: 0.982 sec; CPU: 0.735 sec
[Fri May 14 10:56:57 2021]
Finished job 5.
1 of 7 steps (14%) done

[Fri May 14 10:56:57 2021]
rule sort:
    input: mapped/B.bam
    output: mapped/B.sorted.bam
    jobid: 4
    wildcards: sample=B

Removing temporary output file mapped/B.bam.
[Fri May 14 10:56:58 2021]
Finished job 4.
2 of 7 steps (29%) done

[Fri May 14 10:56:58 2021]
rule bwa:
    input: data/genome.fa, data/samples/A.fastq
    output: mapped/A.bam
    jobid: 3
    wildcards: sample=A

[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 25000 sequences (2525000 bp)...
[M::mem_process_seqs] Processed 25000 reads in 0.707 CPU sec, 0.706 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 1 data/genome.fa data/samples/A.fastq
[main] Real time: 0.962 sec; CPU: 0.732 sec
[Fri May 14 10:57:00 2021]
Finished job 3.
3 of 7 steps (43%) done

[Fri May 14 10:57:00 2021]
rule sort:
    input: mapped/A.bam
    output: mapped/A.sorted.bam
    jobid: 2
    wildcards: sample=A

Removing temporary output file mapped/A.bam.
[Fri May 14 10:57:01 2021]
Finished job 2.
4 of 7 steps (57%) done

[Fri May 14 10:57:01 2021]
rule call:
    input: data/genome.fa, mapped/A.sorted.bam, mapped/B.sorted.bam
    output: calls/all.vcf
    jobid: 1

[warning] samtools mpileup option `g` is functional, but deprecated. Please switch to using bcftools mpileup in future.
[mpileup] 2 samples in 2 input files
Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid
[Fri May 14 10:57:02 2021]
Finished job 1.
5 of 7 steps (71%) done

[Fri May 14 10:57:02 2021]
rule stats:
    input: calls/all.vcf
    output: plots/quals.svg
    jobid: 6

[Fri May 14 10:57:03 2021]
Finished job 6.
6 of 7 steps (86%) done

[Fri May 14 10:57:03 2021]
localrule all:
    input: calls/all.vcf, plots/quals.svg
    jobid: 0

[Fri May 14 10:57:03 2021]
Finished job 0.
7 of 7 steps (100%) done
Complete log: /home/vanessa/Desktop/Code/snakemake/tutorial/.snakemake/log/2021-05-14T105654.746372.snakemake.log
```

After running you will see result data files in newly generated folders mapped, calls, and plots!

```bash
$ tree -L 1
.
├── calls
├── data
├── envs
├── mapped
├── plots
├── README.md
├── scripts
└── Snakefile
```

The resulting plot is kept alongside this repository to show you:

[img/quals.svg](img/quals.svg)

What will happen when you run the workflow with `--use-spack` is that a spack
view will be created for each of your environments. This means that a `spack install`
will be run for each environment to create the executables and other supporting
libraries. We can peek into the hidden .snakemake directory and see our spack
environments:

```bash
$ tree .snakemake/spack
.snakemake/spack
├── 4ff8d68292f72d69cfaa19e8a6ddff03
│   ├── env_setup_done
│   ├── env_setup_start
│   ├── spack.lock
│   └── spack.yaml
├── 6433cc1f09897cb43ae31919ad8e1892
│   ├── env_setup_done
│   ├── env_setup_start
│   ├── spack.lock
│   └── spack.yaml
└── 827d39341264ddaff183ebb496b21eb2
    ├── env_setup_done
    ├── env_setup_start
    ├── spack.lock
    └── spack.yaml
```

Within each folder named by the hash, you see the actual spack view (which
mimics a small filesystem with include, bin, lib, etc.

```bash
.snakemake/spack/827d39341264ddaff183ebb496b21eb2/.spack-env/
├── repos
│   └── builtin
├── transaction_lock
└── view
    ├── bin
    ├── docs
    ├── etc
    ├── include
    ├── lib
    ├── libexec
    ├── man
    ├── sbin
    └── share
```

Note that although spack has support for copy/relocate of environments, the path
is limited in its length, so already in the snakemake metadata directory
with a hash for the name, it's really too long to feasibly try to support.
