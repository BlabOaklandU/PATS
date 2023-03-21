#!/usr/bin/perl

use warnings;
use strict;

my $check_case="";
if(scalar @ARGV == 2){
	$check_case=$ARGV[1];
}

open(FILE,"<$ARGV[0]") || die("Error, could not open file $ARGV[0]: $!");
my @in = <FILE>;
close(FILE);

if($check_case eq "ring4_K5"){
	# check if there are 4 groups with 5 species and 5 proteins (ring of 4*K5)
	exit( 4 != (scalar grep {my @arr=split("\t",$_); scalar @arr >2 && $arr[0] == 5 && $arr[1] == 5 } grep {$_=~m/^[^#]/} @in));
}elsif($check_case eq "P3_K5"){
	# check if there are 3 groups with 5 species and 5 proteins (Path of 3*K5)
	exit( 3 != (scalar grep {my @arr=split("\t",$_); scalar @arr >2 && $arr[0] == 5 && $arr[1] == 5 } grep {$_=~m/^[^#]/} @in));
}elsif($check_case eq "P4"){
	# check if there is 1 groups with all entries
	exit( 1 != (scalar grep {my @arr=split("\t",$_); scalar @arr >2 && $arr[0] == 4 && $arr[1] == 4 } grep {$_=~m/^[^#]/} @in));
}elsif($check_case eq "K5"){
	# check if there is 1 groups with all entries
	exit( 1 != (scalar grep {my @arr=split("\t",$_); scalar @arr >2 && $arr[0] == 5 && $arr[1] == 5 } grep {$_=~m/^[^#]/} @in));
}elsif (scalar(@in) < 15) {
	print STDERR "Outputfile $ARGV[0] is smaller than expected. Test failed!\n";
	exit 1;
}
