#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o /nethome/spusuario14/Entrega2/outputOpenMP.txt
#SBATCH -e /nethome/spusuario14/Entrega2/erroresOpenMP.txt
export OMP_NUM_THREADS=$2
./salidaOpenMp $1 $2
