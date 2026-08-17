#!/usr/bin/env bash
#SBATCH --job-name=deschamps
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#
# Cluster-specific settings (partition, account, qos) are intentionally not
# hard-coded here. Pass them on the command line, e.g.:
#
#   sbatch --partition=<partition> --account=<account> --qos=<qos> deschamps.sh
#
./deschamps
