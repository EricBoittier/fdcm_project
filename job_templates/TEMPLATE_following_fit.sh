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
stepsize={{stepsize}}
learning_rate={{learning_rate}}
skipqs={{skipqs}}
n_charges={{n_charges}}
scan_name={{scan_name}}
suffix={{suffix}}
cubes_dir={{cubes_dir}}
output_dir={{output_dir}}
frames={{frames}}
initial_fit={{initial_fit}}
#initial_fit_cube={{initial_fit_cube}}
prev_frame={{prev_frame}}
start_frame={{start_frame}}
next_frame={{next_frame}}
acd={{acd}}

start=$start_frame
next=$next_frame
dir='frame_'$start"_"$next
output_name=$output_dir/$dir/$dir'-'$start'-'$next'.xyz'
initial_fit=$output_dir/'frame_'$prev_frame'_'$start/'frame_'$prev_frame'_'$start'-'$prev_frame'-'$start'.xyz'
#  Go to the output directory
mkdir -p $output_dir
cd $output_dir

mkdir -p $dir
cd $dir
# Do Initial Fit
#  for initial fit
esp1=$cubes_dir/$scan_name$start$suffix'.p.cube'
dens1=$cubes_dir/$scan_name$start$suffix'.d.cube'
esp=$cubes_dir/$scan_name$next$suffix'.p.cube'
dens=$cubes_dir/$scan_name$next$suffix'.d.cube'

# adjust reference frame
python $ars -charges $initial_fit -pcube $dens1 -pcube2 $dens -frames $frames -output $output_name -acd $acd > $output_name.ARS.log
# do gradient descent fit
$fdcm -xyz $output_name.global -dens $dens -esp  $esp -stepsize $stepsize -n_steps $n_steps -learning_rate $learning_rate -skipqs $skipqs -output $output_name > $output_name.GD.log
# adjust reference frame
python $ars -charges $output_name -pcube $esp  -pcube2 $esp -frames $frames -output $output_name -acd $acd > $output_name.ARS-2.log
# make a cube file for the fit
$cubefit -v -generate -esp $esp -dens $dens  -xyz refined.xyz > $output_name.cubemaking.log
# do analysis
$cubefit -v -analysis -esp $esp -esp2 $n_charges'charges.cube' -dens  $dens > $output_name.analysis.log
echo $PWD





