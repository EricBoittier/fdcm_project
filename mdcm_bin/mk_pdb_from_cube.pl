#!/usr/bin/perl

# Mike Devereux, quick script to read a Gaussian cube file and write a .xyz file

use strict;

if(@ARGV+0 != 1){ die "usage: mk_pdb_from_cube.pl <cube-file> <resname>"; }
my @elements=('H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca');

my $tobohr=1.88973;

my $CUB=$ARGV[0];
my $natm;

# cube file in Gaussian format, needed for the coordinates in the header
open(PDB,">converted.pdb"); # output cube file

print PDB "REMARK       1\nREMARK  pdb generated from $ARGV[0]";

open(REF,"<$CUB");
my $l=0;
my $head=6; # number of lines in header
while(<REF>){

  chomp;
  $l++;
  my $line=$_;
  my @a=split;

  if($l>6 && $l<=$head){
    my $te=$elements[$a[0]-1];
    my $x=$a[2]/$tobohr;
    my $y=$a[3]/$tobohr;
    my $z=$a[4]/$tobohr;
    printf PDB "%s %8.5f %8.5f %8.5f\n",$te,$x,$y,$z;
  }
  if($l==3){
    $head+=$a[0];
    $natm=$a[0];
    print PDB "$natm\n\n";
  }
}
close(REF);

