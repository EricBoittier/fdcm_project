#!/bin/bash
#SBATCH --job-name={{output_dir}}
#SBATCH --nodes=1
#SBATCH --ntasks=1
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
start_frame={{start_frame}}
next_frame={{next_frame}}

start=$start_frame
next=$next_frame
dir='frame_'$start
output_name=$output_dir/$dir/$dir'-'$start'-'$next'.xyz'
#  Go to the output directory
mkdir -p $output_dir
cd $output_dir

mkdir -p $dir
cd $dir
# Do Initial Fit
#  for initial fit
esp=$cubes_dir/$scan_name$next$suffix'.p.cube'
dens=$cubes_dir/$scan_name$next$suffix'.d.cube'
# adjust reference frame
python $ars $initial_fit $initial_fit_cube.d.cube  $esp $frames $output_name > ARS.log
# do gradient descent fit
$fdcm -xyz $output_name.global -dens $dens -esp  $esp -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 -output $output_name > GD.log
python $ars $output_name $esp  $esp $frames $output_name > ARS.log
# make a cube file for the fit
$cubefit -v -generate -esp $esp -dens $dens  -xyz refined.xyz > cubemaking.log
# do analysis
$cubefit -v -analysis -esp $esp -esp2 $n_charges'charges.cube' -dens  $dens > analysis.log