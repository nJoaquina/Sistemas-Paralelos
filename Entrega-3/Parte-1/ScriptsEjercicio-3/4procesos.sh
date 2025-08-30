#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --tasks-per-node=4
#SBATCH -o /nethome/spusuario14/Entrega3/Parte1/outMPI.txt
#SBATCH -e /nethome/spusuario14/Entrega3/Parte1/erroresMPI.txt
mpirun salidaBlocking-ring $1

