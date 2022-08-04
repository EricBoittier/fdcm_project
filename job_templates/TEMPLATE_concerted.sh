#!/bin/bash
#SBATCH --job-name={{output_dir}}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-35
#SBATCH --partition=short
#SBATCH --output={{output_dir}}_%A-%a.out

hostname
#  Path to scripts and executables
cubefit={{cubefit_path}}
fdcm={{fdcm_path}}
ars={{ars_path}}

#  Variables for the job
n_steps={{n_steps}}
n_charges={{n_charges}}
scan_name={{scan_name}}
suffix={{suffix}}
cubes_dir={{cubes_dir}}
output_dir={{output_dir}}
frames={{frames}}
initial_fit={{initial_fit}}
initial_fit_cube={{initial_fit_cube}}
#  for initial fit
esp=$cubes_dir/$scan_name'0'$suffix'.p.cube'
dens=$cubes_dir/$scan_name'0'$suffix'.d.cube'

#  Go to the output directory
mkdir -p $output_dir
cd $output_dir

# Do Initial Fit
#
# adjust reference frame
python $ars $initial_fit $initial_fit_cube.d.cube  $esp $frames 0_fit.xyz > ARS.log
# do gradient descent fit
$fdcm -xyz 0_fit.xyz.global -dens $dens -esp  $esp -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > GD.log
# make a cube file for the fit
$cubefit -v -generate -esp $esp -dens $dens  -xyz refined.xyz > cubemaking.log
# do analysis
$cubefit -v -analysis -esp $esp -esp2 $n_charges'charges.cube' -dens  $dens > analysis.log

initial_fit='../refined.xyz'

#
# Do concerted fit with Slurm array jobs
#
start=$SLURM_ARRAY_TASK_ID
next=$(($start+1))
dir='frame_'$next
output_name=$output_dir/$dir/$dir'-'$start'-'$next'.xyz'

dens=$cubes_dir/$scan_name$start$suffix'.d.cube'
esp=$cubes_dir/$scan_name$start$suffix'.p.cube'

mkdir -p $dir
cd $dir
echo $PWD

python $ars $initial_fit $cubes_dir/$scan_name'0'$suffix'.d.cube' $dens $frames $output_name > ARS.log

cp $output_name'.global' refined.xyz

$fdcm -xyz refined.xyz -dens $dens -esp $esp  -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > GD.log
cp refined.xyz $next'_final.xyz'
# re-adjust to local
python $ars refined.xyz $initial_fit_cube.d.cube $dens $frames refined.xyz > ARS-2.log
# make a cube file for the fit
$cubefit -v -generate -dens $dens -esp $esp  -xyz refined.xyz > cubemaking.log
# do analysis
$cubefit -v -analysis -esp $esp -esp2 $n_charges'charges.cube' -dens  $dens > analysis.log


