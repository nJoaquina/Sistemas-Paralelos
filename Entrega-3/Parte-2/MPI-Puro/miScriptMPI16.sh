#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --tasks-per-node=8
#SBATCH -o /nethome/spusuario14/Entrega3/Parte2/outputMPI.txt
#SBATCH -e /nethome/spusuario14/Entrega3/Parte2/erroresMPI.txt
mpirun salida_mpi $1
