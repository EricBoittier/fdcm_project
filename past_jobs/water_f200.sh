#!/bin/bash
#SBATCH --job-name=/data/unibas/boittier/water_remake200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --output=/home/unibas/boittier/FDCM/out_files/scan_%A.out

hostname
#  Path to scripts and executables
#cubefit=/home/unibas/boittier/fdcm_project/mdcm_bin/cubefit.x
cubefit=/home/unibas/boittier/MDCM/bin/pcubefit.x
fdcm=/home/unibas/boittier/fdcm_project/fdcm.x
ars=/home/unibas/boittier/fdcm_project/ARS.py

#  Variables for the job
n_steps=200
n_charges=6
scan_name=scan
suffix=
cubes_dir=/data/unibas/boittier/models/water/scan
output_dir=/data/unibas/boittier/water_remake200
frames=/home/unibas/boittier/water_model/frames.txt
initial_fit=/home/unibas/boittier/water_model/6-combined_new.xyz
initial_fit_cube=/data/unibas/boittier/models/water/scan/scan0
morton_start=0
acd=/home/unibas/boittier/water_model/water.acd

#  for initial fit
esp=$cubes_dir/$scan_name$morton_start$suffix'.p.cube'
dens=$cubes_dir/$scan_name$morton_start$suffix'.d.cube'

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
$fdcm -xyz refined.xyz.global -dens $dens -esp  $esp -stepsize 0.2 -n_steps 0 -learning_rate 0.5 > GD.log
# re adjust to local
python $ars -charges refined.xyz -pcube $dens -pcube2 $dens -frames $frames -output refined.xyz -acd $acd> ARS-2.log
# make a cubefile for the fit
$cubefit -v -generate -esp $esp -dens $dens  -xyz refined.xyz > cubemaking.log
# do analysis
$cubefit -v -analysis -esp $esp -esp2 $n_charges'charges.cube' -dens  $dens > analysis.log

initial_fit='../frame_'$morton_start/"refined.xyz"
initial_fit_cube=$cubes_dir/$scan_name$morton_start$suffix

cd ..
#
#  Work sequentially through scan
#
start=$morton_start

for next in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
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



