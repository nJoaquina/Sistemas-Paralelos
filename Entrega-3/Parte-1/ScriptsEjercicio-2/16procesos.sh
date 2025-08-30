#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --tasks-per-node=8
#SBATCH -o /nethome/spusuario14/Entrega3/Parte1/outBlocking.txt
#SBATCH -e /nethome/spusuario14/Entrega3/Parte1/erroresBlocking.txt
mpirun salidaBlocking $1
