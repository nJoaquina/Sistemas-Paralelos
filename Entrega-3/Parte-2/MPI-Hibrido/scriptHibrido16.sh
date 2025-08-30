#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive 
#SBATCH --tasks-per-node=1  
#SBATCH -o /nethome/spusuario14/Entrega3/Parte2/hibrido_entrega/outputHibrido.txt 
#SBATCH -e /nethome/spusuario14/Entrega3/Parte2/hibrido_entrega/errorsHibrido.txt 
export OMP_NUM_THREADS=8
mpirun --bind-to none salida_hibrido $1 $2