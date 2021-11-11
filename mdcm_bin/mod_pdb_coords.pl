#!/usr/bin/perl
#
# This script reads a template PDB file with a molecule and dummy "K" atom, and a
# separate xyz file containing new coords for the same molecule without the K atom
# (atom ordering must match), and updates the coords in the template PDB to create
# a new PDB "modified.pdb"

use strict;

if(@ARGV+0 != 2) { die "usage: mod_pdb_coords.pl <pdbfile> <xyzfile>\n"; }

my $pdbfile=$ARGV[0];
my $xyzfile=$ARGV[1];

my @coords;

open(PDB,"<$pdbfile");
open(XYZ,"<$xyzfile");
open(OUT,">$xyzfile.pdb");

print "$pdbfile $xyzfile \n";


my $natm=0;
my $line=0;
while(<XYZ>){

  $line++;
  chomp;
  my @a=split;

  if($line>2 && @a+0 > 3){
    $coords[$natm][0]=$a[1];
    $coords[$natm][1]=$a[2];
    $coords[$natm++][2]=$a[3];
  }

}

close(XYZ);

my $natm2=0;
while(<PDB>){

  chomp;
  my $line=$_;
  my @a=split;

  if($a[0] eq "ATOM" && $a[2] ne "K"){
    printf OUT "ATOM    %3i  %-3s %-4s  %3i      %6.3f  %6.3f  %6.3f  %4.2f  %4.2f           %s\n",
       $a[1],$a[2],$a[3],$a[4],$coords[$natm2][0],$coords[$natm2][1],$coords[$natm2++][2],
       $a[8],$a[9],$a[10];
  }else{
    print OUT "$line\n";
  }

}

if($natm != $natm2){
  print "WARNING: number of atoms in xyz and pdb files seems to differ!\n";
}

