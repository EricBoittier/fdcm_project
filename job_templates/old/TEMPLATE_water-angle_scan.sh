#!/bin/bash
#SBATCH --job-name={{output_dir}}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/home/boittier/FDCM/out_files/{{output_dir}}.out
hostname

olivers_code="/home/boittier/esp_check/bin/cubefit.x"

c_w_dir=$PWD
cd $c_w_dir
cd ..

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
python ~/AdjustReference-System/ARS.py $initial_fit $initial_fit_cube".com.chk.d.cube"  $cubes_dir/$scan_name".com.chk.p.cube" $frames 0_fit.xyz > "ARS.log"
# do gradient descent fit
time ../cubefit.x -xyz global_0_fit.xyz -dens $cubes_dir/$scan_name".com.chk.d.cube" -esp  $cubes_dir/$scan_name".com.chk.p.cube" -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > "GD.log"
# make a cubefile for the fit
$olivers_code -v -generate -esp  $cubes_dir/$scan_name".com.chk.p.cube" -dens  $cubes_dir/$scan_name".com.chk.d.cube" -xyz refined.xyz > "cubemaking.log"
# do analysis
$olivers_code -v -analysis -esp  $cubes_dir/$scan_name".com.chk.p.cube" -esp2 $n_charges"charges.cube" -dens  $cubes_dir/$scan_name".com.chk.d.cube" > "analysis.log"

# do analysis
$olivers_code -v -analysis -esp  $cubes_dir/$scan_name".com.chk.p.cube" -esp2 $initial_fit_cube".com.chk.p.cube" -dens  $cubes_dir/$scan_name".com.chk.d.cube" > "analysis_abinitio.log"


for start in {6..12}
do
next=$(($start+1))
output_name=$output_dir"-"$start"-"$next".xyz"
dir=$start"_"$next
mkdir -p $dir
cd $dir
echo $PWD

python ~/AdjustReference-System/ARS.py $initial_fit $cubes_dir/$start$scan_name"d.cube" $cubes_dir/$next$scan_name"d.cube" $frames $output_name $dih > "ARS.log"

cp "global_"$output_name refined.xyz

time ../../cubefit.x -xyz refined.xyz -dens $cubes_dir/"$next"$scan_name"d.cube" -esp  $cubes_dir/"$next"$scan_name"p.cube" -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > "GD.log"

cp refined.xyz $next"_final.xyz"

# make a cubefile for the fit
$olivers_code -v -generate -esp  $cubes_dir/"$next"$scan_name"p.cube" -dens  $cubes_dir/"$next"$scan_name"d.cube" -xyz refined.xyz > "cubemaking.log"
# do analysis
$olivers_code -v -analysis -esp  $cubes_dir/"$next"$scan_name"p.cube" -esp2 $n_charges"charges.cube" -dens  $cubes_dir/"$next"$scan_name"d.cube" > "analysis.log"

initial_fit=../$start"_"$next"/"$next"_final.xyz"
start=$(($start+1))
next=$(($start+1))

cd ..
done


