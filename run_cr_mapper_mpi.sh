#!/bin/bash
#PBS -d ${WORKDIR}
#PBS -l nodes=${NODEs}:ppn=${PPN}
#PBS -N cr_mapper
#PBS -l walltime=${WALLTIME}
#PBS -l pmem=${MEM}
#PBS -q ${QUEUE}

# This and the companion python script implements the cr mapper over a list
# ${IN_FILE_LIST} of reference source files across multiple compute nodes using
# mpi and sends solution files to ${OUT_FILE_LIST}.

# cr_mapper.py when tested used memory at about 10x the size of the reference
# source file, but test this on the deployment system.

# Launch using the hosts
mpirun --hostfile $PBS_NODEFILE python3 -m mpi4py.futures run_cr_mapper_mpi.py ${IN_FILE_LIST} ${OUT_FILE_LIST}
