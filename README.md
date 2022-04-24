# SnapperRocks

Analysis pipeline by Forde Lab.

Bioinfomatics development by Dr. Brian Forde, coding by Thom Cuddihy

## Requirements

SnapperRocks requires that [Nextflow](https://www.nextflow.io/) v20+ be installed to run, which is compatible with Linux and OSX.

Compute environments are provided in [Singularity/Apptainer](https://apptainer.org/) containers for local and cluster deployment. For AWS deployment, instances may be created from the containers.

## Containers

After pulling, you will need to download the latest published Singularity containers:

* [General Container](https://fordelab.com/share/sr_general_db_0.9.img), MD5: eb324ed5ee10d303a2b2d92358565388
* [Kraken2 Container](https://fordelab.com/share/sr_kraken_db_0.9.img), MD5: fac635cecaf4af8f28c8b01252b19900
* [Nesoni Container](https://fordelab.com/share/sr_nesoni_0.9.img), MD5: df1b5c1deef639ccdae05a2d66364ed1

## Initial Configuration

Update the `nextflow.config` file as appropriate for Singularity container location (default is '/data/images/').

## Running

SnapperRocks may be launched from the command line, using the Nextflow binary and the pipeline’s main.nf file. Any changes to the default values (see User Guide) may be specified with a `--` prefix and then the value. For Boolean parameters, a true value may be indicated with just the `--` prefix. 

Defaults may be permanently changed by editing the nextflow.config file in the path of SnapperRocks’ main.nf.

E.g.:
```bash
nextflow /path/to/SnapperRocks/main.nf --run_id "run_71" --fastq "/ftp/ingest/*_R{1,2}.fastq.gz" --trim --taxoprofile  --assembly --typing --nesoni --cluster --executor "slurm"
```

At a minimum, the values that change between runs should be specified on the command line (e.g. `run_id`, and depending on data ingest and handling `fastq` and/or `results`)

Parameters that will stay the same between runs may be updated in the nextflow.config or left as command line arguments and stored as a batch file for record keeping.

## Example Data

Example data is available for download [here](https://fordelab.com/share/Example_Data.zip), MD5: 556e88b1b454d40d17e24376da2e1e43