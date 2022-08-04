#!/bin/bash
#SBATCH --job-name=ester
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/home/boittier/ester_traj.out
hostname

olivers_code="/home/boittier/esp_check/bin/cubefit.x"
pcubex="/home/boittier/FDCM/cubefit.x"

c_w_dir=$PWD
cd $c_w_dir
cd ..

n_steps=10
n_charges=24
scan_name="frame_"
cubes_dir="/home/boittier/ester_traj/t1"
output_dir="ester_t2_"$n_steps
frames="/home/boittier/ester_traj/frames.txt"
initial_fit="/home/boittier/FDCM/models/ester/0_1/refined.xyz"
initial_fit_cube="/home/boittier/cube_files/ester/1_SCAN_1_2_6_8_F/O"

starting_frame=0

first_fit=$cubes_dir/"frame_"$starting_frame.chk

mkdir -p $output_dir
cd $output_dir

# Do Initial Fit to 
mkdir "frame_"$starting_frame
cd "frame_"$starting_frame

# adjust reference frame
echo $initial_fit $initial_fit_cube".d.cube"  $first_fit.p.cube $frames first_fit.xyz

python ~/AdjustReference-System/ARS.py $initial_fit $initial_fit_cube".d.cube"  $first_fit.p.cube $frames fit.xyz > "ARS1.log"

# do gradient descent fit
time $pcubex -xyz global_fit.xyz -dens $first_fit.d.cube -esp  $first_fit.p.cube -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > "GD.log"

python ~/AdjustReference-System/ARS.py refined.xyz $first_fit.p.cube $first_fit.p.cube $frames 6_fit.xyz > "ARS.log"

# make a cubefile for the fit
$olivers_code -v -generate -esp  $first_fit.p.cube -dens  $first_fit.d.cube -xyz refined.xyz > "cubemaking.log"

# do analysis
$olivers_code -v -analysis -esp  $first_fit.p.cube -esp2 $n_charges"charges.cube" -dens $first_fit.d.cube > "analysis.log"

cd ..

initial_fit="../frame_"$starting_frame"/refined.xyz"

for start in {0..498}
do
start=$(($start + $starting_frame))
next=$(($start + 1))
output_name=$output_dir"-"$start"-"$next".xyz"
dir="frame_"$next

startd="frame_"$start".chk.d.cube"
startp="frame_"$start".chk.p.cube"
nextd="frame_"$next".chk.d.cube"
nextp="frame_"$next".chk.p.cube"

mkdir -p $dir
cd $dir
echo $PWD

python ~/AdjustReference-System/ARS.py $initial_fit $cubes_dir/$startd $cubes_dir/$nextd $frames $output_name $dih > "ARS.log"

cp "global_"$output_name refined.xyz

time $pcubex -xyz refined.xyz -dens $cubes_dir/$nextd -esp  $cubes_dir/$nextp -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > "GD.log"

python ~/AdjustReference-System/ARS.py refined.xyz $cubes_dir/$nextd $cubes_dir/$nextd $frames $output_name $dih > "ARS.log"

cp refined.xyz $next"_final.xyz"

# make a cubefile for the fit
$olivers_code -v -generate -esp  $cubes_dir/$nextp  -dens $cubes_dir/$nextd  -xyz refined.xyz > "cubemaking.log"
# do analysis
$olivers_code -v -analysis -esp  $cubes_dir/$nextp -esp2 $n_charges"charges.cube" -dens  $cubes_dir/$nextd > "analysis.log"

initial_fit=../$dir"/"$next"_final.xyz"
start=$(($start+1))
next=$(($start+1))

cd ..
done


