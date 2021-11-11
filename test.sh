#!/bin/bash
#SBATCH --job-name=/home/boittier/fdcm_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --output=/home/unibas/boittier/FDCM/out_files/SCAN_amide1.pdb-_%A-%a.out

hostname
#  Path to scripts and executables
cubefit=/home/boittier/Documents/PhD/fdcm_project/mdcm_bin/cubefit.x
fdcm=/home/boittier/Documents/PhD/fdcm_project/fdcm.x
ars=/home/boittier/Documents/PhD/fdcm_project/ARS.py

#  Variables for the job
n_steps=1
n_charges=20
scan_name=SCAN_amide1.pdb-
suffix=.xyz.chk
cubes_dir=/home/boittier/Documents/data
output_dir=/home/boittier/fdcm_test
frames=/home/boittier/Documents/PhD/fdcm_project/mdcms/amide/model2/frames.txt
initial_fit=/home/boittier/Documents/PhD/fdcm_project/mdcms/amide/model2/20_charges_refined.xyz
initial_fit_cube=/home/boittier/Documents/PhD/fdcm_project/mdcms/amide/model2/1_1_new-0.chk
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
# make a cubefile for the fit
$cubefit -v -generate -esp $esp -dens $dens  -xyz refined.xyz > cubemaking.log
# do analysis
$cubefit -v -analysis -esp $esp -esp2 $n_charges'charges.cube' -dens  $dens > analysis.log

initial_fit='../refined.xyz'

#
#  Work sequentially through scan
#
for start in {0..10}
do
start=$(($start))
next=$(($start+1))
dir='frame_'$next
output_name=$output_dir/$dir/$dir'-'$start'-'$next'.xyz'

dens=$cubes_dir/$scan_name$start$suffix'.d.cube'
esp=$cubes_dir/$scan_name$start$suffix'.p.cube'

mkdir -p $dir
cd $dir
echo $PWD

python $ars $initial_fit $initial_fit_cube.d.cube $dens $frames $output_name > ARS.log
cp $output_name'.global' refined.xyz
$fdcm -xyz refined.xyz -dens $dens -esp $esp  -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > GD.log
cp refined.xyz $next'_final.xyz'
# re-adjust to local
python $ars refined.xyz $initial_fit_cube.d.cube $dens $frames refined.xyz > ARS-2.log
# make a cubefile for the fit
$cubefit -v -generate -dens $dens -esp $esp  -xyz refined.xyz > cubemaking.log
# do analysis
$cubefit -v -analysis -esp $esp -esp2 $n_charges'charges.cube' -dens  $dens > analysis.log

initial_fit=../$dir"/"$next"_final.xyz"

cd ..
done



