#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o /nethome/spusuario14/Entrega2/outputSecuencial.txt
#SBATCH -e /nethome/spusuario14/Entrega2/erroresSecuencial.txt
./salidaSecuencial $1
