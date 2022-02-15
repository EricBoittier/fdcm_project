#!/bin/bash
#SBATCH --job-name={{output_dir}}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --output=/home/unibas/boittier/FDCM/out_files/{{scan_name}}_%A.out

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
morton_start={{morton_start}}
acd={{acd}}
#  for initial fit
esp=$cubes_dir/$scan_name'0'$suffix'.p.cube'
dens=$cubes_dir/$scan_name'0'$suffix'.d.cube'

#  Go to the output directory
mkdir -p $output_dir
cd $output_dir || return

mkdir 'frame_'$morton_start
cd 'frame_'$morton_start || return

# Do Initial Fit
#
# adjust reference frame
python $ars -charges $initial_fit -pcube $initial_fit_cube.d.cube  -pcube2 $esp -frames $frames -output refined.xyz -acd $acd> ARS.log
# do gradient descent fit
$fdcm -xyz refined.xyz.global -dens $dens -esp  $esp -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > GD.log
# re adjust to local
python $ars -charges refined.xyz -pcube $dens -pcube2 $dens -frames $frames -output refined.xyz -acd $acd> ARS-2.log
# make a cubefile for the fit
$cubefit -v -generate -esp $esp -dens $dens  -xyz refined.xyz > cubemaking.log
# do analysis
$cubefit -v -analysis -esp $esp -esp2 $n_charges'charges.cube' -dens  $dens > analysis.log

initial_fit='../frame_'$morton_start/"refined.xyz"
initial_fit_cube=$cubes_dir/$scan_name'0'$suffix

cd ..
#
#  Work sequentially through scan
#
start=$morton_start

for next in {{ morton }}
do

dir='frame_'$next
output_name=$output_dir/$dir/$dir'-'$start'-'$next'.xyz'

dens=$cubes_dir/$scan_name$next$suffix'.d.cube'
esp=$cubes_dir/$scan_name$next$suffix'.p.cube'

mkdir -p $dir
cd $dir || return
echo $PWD

python $ars -charges $initial_fit -pcube $initial_fit_cube.d.cube -pcube2 $dens -frames $frames -output $output_name -acd $acd> ARS.log

cp $output_name'.global' refined.xyz
$fdcm -xyz refined.xyz -dens $dens -esp $esp  -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > GD.log
cp refined.xyz $next'_final.xyz'

# re-adjust to local
python $ars -charges refined.xyz -pcube $initial_fit_cube.d.cube -pcube2 $dens -frames $frames -charges refined.xyz -acd $acd> ARS-2.log
# make a cubefile for the fit
$cubefit -v -generate -dens $dens -esp $esp  -xyz refined.xyz > cubemaking.log
# do analysis
$cubefit -v -analysis -esp $esp -esp2 $n_charges'charges.cube' -dens  $dens > analysis.log

initial_fit=../$dir"/"$next"_final.xyz"
initial_fit_cube=$cubes_dir/$scan_name$next$suffix

start=$next

cd ..
done



