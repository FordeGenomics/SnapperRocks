# SnapperRocks

Analysis pipeline by Forde Lab.

Bioinfomatics development by Dr. Brian Forde, coding by Thom Cuddihy


## Installation

After pulling, you will need to merge the taxonomic database back together using `merge_db.sh`

## Misc

After adding new MLST refs, make sure that the bowtie2 db is built for each.

```bash
cd db/MLST

# find refs w/o bt2 db
for d in `ls -d */`; do 
    S="${d/\//}"
    ls ${S}/*.bt2 > /dev/null; 
done

# example bulk generation
for d in `ls -d */`; do
  S="${d/\//}"
  bowtie2-build "/opt/SnapperRocks/db/MLST/${S}/${S}.fa" "/opt/SnapperRocks/db/MLST/${S}/${S}.fa"
done
```

## Docker

We use .dockerignore trickery when building the docker containers. Set working directory to `./docker` and run the build script (with version tag) from there. E.g.:

```bash
cd docker
./nesoni/build.sh 0.8
```

## Singularity

We support singularity containers by building the docker containers first, exporting them to tar using `docker save` and then re-building into Singularity with `singularity build container.img docker-archive://container.tar`

## NextFlow

Update the `nextflow.config` file as appropriate for local system (docker/singularity/native, slurm/pbs/local/awsbatch e.g.) and run options:

```bash
cd NextFlow
vi nextflow.config
nextflow run main.nf
```