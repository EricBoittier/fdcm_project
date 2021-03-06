#!/bin/bash
#SBATCH --job-name={{output_dir}}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-35
#SBATCH --partition=short
#SBATCH --output=/home/unibas/boittier/FDCM/out_files/{{output_dir}}_%A-%a.out

hostname
cubefit="{{cubefit_path}}"
fdcm="{{fdcm_path}}"
ars="{{ars_path}}"
c_w_dir=$PWD


dih="{{dih}}"
n_steps={{n_steps}}
n_charges={{n_charges}}
scan_name="{{scan_name}}"
cubes_dir="{{cubes_dir}}"
output_dir="{{output_dir}}"
frames="{{frames}}"
initial_fit="{{initial_fit}}"
initial_fit_cube="{{initial_fit_cube}}"


mkdir -p $output_dir
cd $output_dir

# Do Initial Fit
#
# adjust reference frame
python $ars $initial_fit $initial_fit_cube.d.cube  $cubes_dir/0$scan_name"p.cube" $frames 0_fit.xyz $dih > "ARS.log"
# do gradient descent fit
$fdcm -xyz 0_fit.xyz.global -dens $cubes_dir/0"$scan_name"d.cube -esp  $cubes_dir/"0"$scan_name"p.cube" -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > "GD.log"
# make a cubefile for the fit
$cubefit -v -generate -esp  $cubes_dir/"0"$scan_name"p.cube" -dens  $cubes_dir/"0"$scan_name"d.cube" -xyz refined.xyz > "cubemaking.log"
# do analysis
$cubefit -v -analysis -esp  $cubes_dir/"0"$scan_name"p.cube" -esp2 $n_charges"charges.cube" -dens  $cubes_dir/"0"$scan_name"d.cube" > "analysis.log"

initial_fit="../refined.xyz"

# 
# Do concerted fit with Slurm array jobs
#
start=$SLURM_ARRAY_TASK_ID
next=$(($start+1))
output_name=$output_dir"-"$start"-"$next".xyz"
dir=$start"_"$next
mkdir -p $dir
cd $dir
echo $PWD

python ~/AdjustReference-System/ARS.py $initial_fit $initial_fit_cube.d.cube $cubes_dir/$start$scan_name"d.cube" $frames $output_name $dih > "ARS.log"
cp $output_name".global" refined.xyz
time ../../cubefit.x -xyz refined.xyz -dens $cubes_dir/"$start"$scan_name"d.cube" -esp  $cubes_dir/"$start"$scan_name"p.cube" -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > "GD.log"

cp refined.xyz $next"_final.xyz"

# re-adjust to local
python ~/AdjustReference-System/ARS.py refined.xyz $initial_fit_cube.d.cube $cubes_dir/$start$scan_name"d.cube" $frames refined.xyz $dih > "ARS-2.log"
# make a cubefile for the fit
$cubefit -v -generate -esp  $cubes_dir/"$start"$scan_name"p.cube" -dens  $cubes_dir/"$start"$scan_name"d.cube" -xyz refined.xyz > "cubemaking.log"
# do analysis
$cubefit -v -analysis -esp  $cubes_dir/"$start"$scan_name"p.cube" -esp2 $n_charges"charges.cube" -dens  $cubes_dir/"$start"$scan_name"d.cube" > "analysis.log"


