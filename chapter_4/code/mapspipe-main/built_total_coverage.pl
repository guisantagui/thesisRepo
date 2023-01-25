#!/usr/bin/env perl 
use strict;use warnings;

my $position=1;
my $end=4411532;

open FILE, $ARGV[0] or die "Can't open Input: $!\n";


while (<FILE>) {
chomp;
my @fields= split ( /\t/, $_);
LI: if(($fields[0]>$position)&&($position<$end)){
	print "$position\t0\n";
	$position++;
	goto LI; 
}

else {
	$position++;	
	print "$_\n";

}

}
while ($position<=$end){
	print "$position\t0\n";
	$position++;
}
exit;
