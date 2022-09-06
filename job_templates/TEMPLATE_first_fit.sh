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
acd={{acd}}
#  Go to the output directory
mkdir -p $output_dir
cd $output_dir

start=$start_frame
next=$next_frame
dir='frame_'$start'_'$next
mkdir -p $dir
cd $dir
# Do Initial Fit
#  for initial fit
output_name=$output_dir/$dir/$dir'-'$start'-'$next'.xyz'
esp=$cubes_dir/$scan_name$start$suffix'.p.cube'
dens=$cubes_dir/$scan_name$start$suffix'.d.cube'
# adjust reference frame
python $ars -charges $initial_fit -pcube $initial_fit_cube.d.cube  -pcube2 $esp -frames $frames -output $start"_fit.xyz" -acd $acd > $output_name.ARS.log
# do gradient descent fit
$fdcm -xyz $start"_fit.xyz.global" -dens $dens -esp  $esp -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 -output $output_name  > $output_name.GD.log
# re-adjust to local
python $ars -charges $output_name -pcube $initial_fit_cube.d.cube -pcube2 $dens -frames $frames -output $output_name -acd $acd > $output_name.ARS-2.log
# make a cube file for the fit
$cubefit -v -generate -esp $esp -dens $dens  -xyz refined.xyz > $output_name.cubemaking.log
# do analysis
$cubefit -v -analysis -esp $esp -esp2 $n_charges'charges.cube' -dens  $dens > $output_name.analysis.log

initial_fit=$output_name
cd ..

dir='frame_'$start'_'$next
output_name=$output_dir/$dir/$dir'-'$start'-'$next'.xyz'
dens=$cubes_dir/$scan_name$next$suffix'.d.cube'
esp=$cubes_dir/$scan_name$next$suffix'.p.cube'

mkdir -p $dir
cd $dir

# Adjust reference frame
python $ars -charges $initial_fit -pcube $cubes_dir/$scan_name$start$suffix'.d.cube' -pcube2 $dens -frames $frames -output $output_name -acd $acd > $output_name.ARS.log
cp $output_name'.global' refined.xyz
$fdcm -xyz refined.xyz -dens $dens -esp $esp  -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 -output $output_name > $output_name.GD.log
cp refined.xyz $next'_final.xyz'
# re-adjust to local
python $ars -charges $output_name -pcube $dens -pcube2 $dens -frames $frames -output $output_name -acd $acd > $output_name.ARS-2.log
# make a cube file for the fit
$cubefit -v -generate -dens $dens -esp $esp  -xyz refined.xyz > $output_name.cubemaking.log
# do analysis
$cubefit -v -analysis -esp $esp -esp2 $n_charges'charges.cube' -dens  $dens > $output_name.analysis.log

echo $PWD




