#!/bin/bash
#SBATCH --job-name=water_angle
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/home/boittier/FDCM/out_files/{{output_dir}}.out
hostname

olivers_code="/home/boittier/esp_check/bin/cubefit.x"

c_w_dir=$PWD
cd $c_w_dir
cd ..

n_steps=0
n_charges=6
scan_name="angle_"
cubes_dir="/home/boittier/v1"
output_dir="angle"
frames="/home/boittier/v1/frames.txt"
initial_fit="/home/boittier/v1/6_charges_refined.xyz"
initial_fit_cube="/home/boittier/v1/2_2_2_new"


mkdir -p $output_dir
cd $output_dir



# Do Initial Fit to 
mkdir "angle_6"
cd "angle_6"

# adjust reference frame
echo $initial_fit $initial_fit_cube".com.chk.d.cube"  $cubes_dir/"angle_6_new.com.chk.fchk.p.cube" $frames 6_fit.xyz

python ~/AdjustReference-System/ARS.py $initial_fit $initial_fit_cube".com.chk.d.cube"  $cubes_dir/"angle_6_new.com.chk.fchk.p.cube" $frames 6_fit.xyz > "ARS.log"

# do gradient descent fit
time ../cubefit.x -xyz global_6_fit.xyz -dens $cubes_dir/"angle_6_new.com.chk.fchk.d.cube" -esp  $cubes_dir/"angle_6_new.com.chk.fchk.p.cube" -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > "GD.log"

python ~/AdjustReference-System/ARS.py refined.xyz $cubes_dir/"angle_6_new.com.chk.fchk.p.cube"  $cubes_dir/"angle_6_new.com.chk.fchk.p.cube" $frames 6_fit.xyz > "ARS.log"

# make a cubefile for the fit
$olivers_code -v -generate -esp  $cubes_dir/"angle_6_new.com.chk.fchk.p.cube" -dens  $cubes_dir/"angle_6_new.com.chk.fchk.d.cube" -xyz refined.xyz > "cubemaking.log"

# do analysis
$olivers_code -v -analysis -esp  $cubes_dir/"angle_6_new.com.chk.fchk.p.cube" -esp2 $n_charges"charges.cube" -dens $cubes_dir/"angle_6_new.com.chk.fchk.d.cube" > "analysis.log"

cd ..

# Back
for start in {6..0}
do

next=$(($start-1))

output_name=$output_dir"-"$start"-"$next".xyz"

dir="angle_"$next

startd="angle_"$start"_new.com.chk.fchk.d.cube"
startp="angle_"$start"_new.com.chk.fchk.p.cube"
nextd="angle_"$next"_new.com.chk.fchk.d.cube"
nextp="angle_"$next"_new.com.chk.fchk.p.cube"

mkdir -p $dir
cd $dir
echo $PWD

initial_fit="../refined.xyz"

python ~/AdjustReference-System/ARS.py $initial_fit $cubes_dir/$startd $cubes_dir/$nextd $frames $output_name $dih > "ARS.log"

cp "global_"$output_name refined.xyz

time ../../cubefit.x -xyz refined.xyz -dens $cubes_dir/$nextd -esp  $cubes_dir/$nextp -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > "GD.log"

cp refined.xyz $next"_final.xyz"

python ~/AdjustReference-System/ARS.py refined.xyz $cubes_dir/$nextd $cubes_dir/$nextd $frames $output_name $dih > "ARS.log"


# make a cubefile for the fit
$olivers_code -v -generate -esp  $cubes_dir/$nextp  -dens $cubes_dir/$nextd  -xyz refined.xyz > "cubemaking.log"
# do analysis
$olivers_code -v -analysis -esp  $cubes_dir/$nextp -esp2 $n_charges"charges.cube" -dens  $cubes_dir/$nextd > "analysis.log"

initial_fit=../$dir"/"$next"_final.xyz"
start=$(($start+1))
next=$(($start+1))

cd ..
done

#Forward

for start in {6..12}
do
next=$(($start+1))
output_name=$output_dir"-"$start"-"$next".xyz"
dir="angle_"$next

startd="angle_"$start"_new.com.chk.fchk.d.cube"
startp="angle_"$start"_new.com.chk.fchk.p.cube"
nextd="angle_"$next"_new.com.chk.fchk.d.cube"
nextp="angle_"$next"_new.com.chk.fchk.p.cube"

mkdir -p $dir
cd $dir
echo $PWD

python ~/AdjustReference-System/ARS.py $initial_fit $cubes_dir/$startd $cubes_dir/$nextd $frames $output_name $dih > "ARS.log"

cp "global_"$output_name refined.xyz

time ../../cubefit.x -xyz refined.xyz -dens $cubes_dir/$nextd -esp  $cubes_dir/$nextp -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > "GD.log"

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


