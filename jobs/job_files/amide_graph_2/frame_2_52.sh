#!/bin/bash
#SBATCH --job-name=/data/unibas/boittier/amide_graph_2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --output=/data/unibas/boittier/amide_graph_2_%A-%a.out

hostname
#  Path to scripts and executables
cubefit=/home/unibas/boittier/fdcm_project/mdcm_bin/cubefit.x
fdcm=/home/unibas/boittier/fdcm_project/fdcm.x
ars=/home/unibas/boittier/fdcm_project/ARS.py
#  Variables for the job
n_steps=2
n_charges=24
scan_name=frame_
suffix=.chk
cubes_dir=/data/unibas/boittier/fdcm/amide_graph
output_dir=/data/unibas/boittier/amide_graph_2
frames=/home/unibas/boittier/fdcm_project/mdcms/amide/model1/frames.txt
initial_fit=/home/unibas/boittier/fdcm_project/mdcms/amide/model1/24_charges_refined.xyz
initial_fit_cube=/home/unibas/boittier/fdcm_project/mdcms/amide/model1/amide1.pdb.chk
start_frame=2
next_frame=52
acd=/home/unibas/boittier/fdcm_project/0_fit.xyz.acd
#  Go to the output directory
mkdir -p $output_dir
cd $output_dir

start=$start_frame
next=$next_frame
dir='frame_'$start
mkdir -p $dir
cd $dir
# Do Initial Fit
#  for initial fit
output_name=$output_dir/$dir/$dir'-'$start'-'$next'.xyz'
esp=$cubes_dir/$scan_name'0'$suffix'.p.cube'
dens=$cubes_dir/$scan_name'0'$suffix'.d.cube'
# adjust reference frame
python $ars -charges $initial_fit -pcube $initial_fit_cube.d.cube  -pcube2 $esp -frames $frames -output 0_fit.xyz -acd $acd > $output_name.ARS.log
# do gradient descent fit
$fdcm -xyz 0_fit.xyz.global -dens $dens -esp  $esp -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 -output $output_name  > $output_name.GD.log
# re-adjust to local
python $ars -charges $output_name -pcube $initial_fit_cube.d.cube -pcube2 $dens -frames $frames -output $output_name -acd $acd > $output_name.ARS-2.log
# make a cube file for the fit
$cubefit -v -generate -esp $esp -dens $dens  -xyz refined.xyz > $output_name.cubemaking.log
# do analysis
$cubefit -v -analysis -esp $esp -esp2 $n_charges'charges.cube' -dens  $dens > $output_name.analysis.log

initial_fit=$output_name
cd ..

dir='frame_'$next
output_name=$output_dir/$dir/$dir'-'$start'-'$next'.xyz'
dens=$cubes_dir/$scan_name$next$suffix'.d.cube'
esp=$cubes_dir/$scan_name$next$suffix'.p.cube'

mkdir -p $dir
cd $dir

# Adjust reference frame
python $ars -charges $initial_fit -pcube $cubes_dir/$scan_name'0'$suffix'.d.cube' -pcube2 $dens -frames $frames -output $output_name -acd $acd > $output_name.ARS.log
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




sbatch /home/unibas/boittier/fdcm_project/job_files/amide_graph_2/frame_2_0.sh 

sbatch /home/unibas/boittier/fdcm_project/job_files/amide_graph_2/frame_2_9.sh 

sbatch /home/unibas/boittier/fdcm_project/job_files/amide_graph_2/frame_2_31.sh 

sbatch /home/unibas/boittier/fdcm_project/job_files/amide_graph_2/frame_2_42.sh 

sbatch /home/unibas/boittier/fdcm_project/job_files/amide_graph_2/frame_2_59.sh 

sbatch /home/unibas/boittier/fdcm_project/job_files/amide_graph_2/frame_2_64.sh 

sbatch /home/unibas/boittier/fdcm_project/job_files/amide_graph_2/frame_52_61.sh 

sbatch /home/unibas/boittier/fdcm_project/job_files/amide_graph_2/frame_52_63.sh 

sbatch /home/unibas/boittier/fdcm_project/job_files/amide_graph_2/frame_52_45.sh 
