#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o /nethome/spusuario14/Entrega2/outputPThreads.txt
#SBATCH -e /nethome/spusuario14/Entrega2/erroresPThreads.txt
./salidaPthread $1 $2
