#!/bin/bash
source /app/lmod/lmod/init/profile
module purge
module load Apptainer/1.1.6

echo 'Running apptainer' >> script_using_mclapply.Rrender.log.txt

apptainer run \
    --bind /fh/fast:/fh/fast \
    https://sif-registry.fredhutch.org/bioconductor_docker_RELEASE_3_22-R-4.5.2.sif \
    bash script_using_mclapply.Rrender.Rcode.sh >> script_using_mclapply.Rrender.log.txt

echo 'Finished running apptainer' >> script_using_mclapply.Rrender.log.txt


module purge
