#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $sizes = shift;
my $dat = shift;

die "No sizes data" unless $dat && -s $dat;

my %hash;

open (IN,$sizes);
while (my $ln=<IN>){
	chomp($ln);
	next unless $ln;
	my @d = split("\t",$ln);
	$hash{$d[0]} = $d[1];
}

close IN;
open (IN2,$dat);
while (my $ln=<IN2>){
        chomp($ln);
	next unless $ln;
	if ($ln=~/^@(\S+)/){
		my $id = $1;
		if ($ln=~/fractionMasked=([\d\.]+)/){
			my $fraction = $1;
			print $id."\t".$hash{$id}."\t".$fraction."\n";
		}
	}
}
