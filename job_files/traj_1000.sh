#!/bin/bash
#SBATCH --job-name=amide
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/home/boittier/FDCM/out_files/amide.out
hostname

olivers_code="/home/boittier/esp_check/bin/cubefit.x"
pcubex="/home/boittier/FDCM/cubefit.x"

c_w_dir=$PWD
cd $c_w_dir
cd ..

n_steps=0
n_charges=20
scan_name="frame_"
cubes_dir="/home/boittier/amide_dih/scan"
output_dir="amide_dih_0"
frames="/data/boittier/pydcm/amide1/frames.txt"
initial_fit="/home/boittier/pydcm/amide1/5-refine/20_charges_refined.xyz"
initial_fit_cube="/home/boittier/amide_dih/OPT_B3LYP_6-311++G**/amide1.pdb.chk"

starting_frame=499
first_fit=$cubes_dir/"frame_"$starting_frame"/frame_"$starting_frame

mkdir -p $output_dir
cd $output_dir

# Do Initial Fit to 
mkdir "frame_"$starting_frame
cd "frame_"$starting_frame

# adjust reference frame
echo $initial_fit $initial_fit_cube".d.cube"  $first_fit.p.cube $frames first_fit.xyz

python ~/AdjustReference-System/ARS.py $initial_fit $initial_fit_cube".d.cube"  $first_fit.p.cube $frames fit.xyz > "ARS.log"

# do gradient descent fit
time $pcubex -xyz global_fit.xyz -dens $first_fit.d.cube -esp  $first_fit.p.cube -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > "GD.log"

python ~/AdjustReference-System/ARS.py refined.xyz $first_fit.p.cube $first_fit.p.cube $frames 6_fit.xyz > "ARS.log"

# make a cubefile for the fit
$olivers_code -v -generate -esp  $first_fit.p.cube -dens  $first_fit.d.cube -xyz refined.xyz > "cubemaking.log"

# do analysis
$olivers_code -v -analysis -esp  $first_fit.p.cube -esp2 $n_charges"charges.cube" -dens $first_fit.d.cube > "analysis.log"

cd ..

initial_fit="../frame_"$starting_frame"/refined.xyz"

for start in {0..100}
do
start=$(($start + $starting_frame))
next=$(($start + 1))
output_name=$output_dir"-"$start"-"$next".xyz"
dir="frame_"$next

startd="frame_"$start"/frame_"$start".d.cube"
startp="frame_"$start"/frame_"$start".p.cube"
nextd="frame_"$next"/frame_"$next".d.cube"
nextp="frame_"$next"/frame_"$next".p.cube"

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


