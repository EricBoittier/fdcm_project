#!/bin/bash
#SBATCH --job-name=/data/unibas/boittier/morton-large-control5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --output=/home/unibas/boittier/FDCM/out_files/SCAN_amide1.pdb-_%A-%a.out

hostname
#  Path to scripts and executables
cubefit=/home/unibas/boittier/fdcm_project/mdcm_bin/cubefit.x
fdcm=/home/unibas/boittier/fdcm_project/fdcm.x
ars=/home/unibas/boittier/fdcm_project/ARS.py

#  Variables for the job
n_steps=5
n_charges=24
scan_name=SCAN_amide1.pdb-
suffix=.xyz.chk
cubes_dir=/data/unibas/boittier/fdcm/amide/scan-large
output_dir=/data/unibas/boittier/morton-large-control5
frames=/home/unibas/boittier/fdcm_project/mdcms/amide/model1/frames.txt
initial_fit=/home/unibas/boittier/fdcm_project/mdcms/amide/model1/24_charges_refined.xyz
initial_fit_cube=/home/unibas/boittier/fdcm_project/mdcms/amide/model1/amide1.pdb.chk
morton_start=0
#  for initial fit
esp=$cubes_dir/$scan_name'0'$suffix'.p.cube'
dens=$cubes_dir/$scan_name'0'$suffix'.d.cube'

#  Go to the output directory
mkdir -p $output_dir
cd $output_dir

mkdir 'frame_'$morton_start
cd 'frame_'$morton_start

# Do Initial Fit
#
# adjust reference frame
python $ars $initial_fit $initial_fit_cube.d.cube  $esp $frames refined.xyz > ARS.log
# do gradient descent fit
$fdcm -xyz refined.xyz.global -dens $dens -esp  $esp -stepsize 0.2 -n_steps $n_steps -learning_rate 0.5 > GD.log
# re adjust to local
python $ars refined.xyz $dens $dens $frames refined.xyz > ARS-2.log
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

for next in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159
do

dir='frame_'$next
output_name=$output_dir/$dir/$dir'-'$start'-'$next'.xyz'

dens=$cubes_dir/$scan_name$next$suffix'.d.cube'
esp=$cubes_dir/$scan_name$next$suffix'.p.cube'

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
initial_fit_cube=$cubes_dir/$scan_name$next$suffix

start=$next

cd ..
done



