#!/usr/bin/env bash
#SBATCH --job-name=deschamps
#SBATCH --partition=etna-shared
#SBATCH --account=nano
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
./deschamps
