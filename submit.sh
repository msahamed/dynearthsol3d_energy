#!/bin/sh
#PBS -l nodes=1:default:ppn=1
#PBS -l walltime=600:00:00
#PBS -A [CERI]
#PBS -N diffusion_combined

#export OMP_NUM_THREADS=1

# source the module command

source /etc/profile.d/modules.sh

cd /home/msahamed/nmdExperiemnt_combined/diffusion_DES/result

declare -a arr=("core_complex_deviatoric_variable_dilation" "core_complex_energy_fixed_dilation" "core_complex_energy_variable_dilation" "core_complex_energy")
for i in "${arr[@]}"
do
    file_name="$i.sh"
    qsub file_name
done
