#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o /nethome/spusuario14/output.txt
#SBATCH -e /nethome/spusuario14/errors.txt
./salidaPunto2 $1 $2
