#!/usr/bin/perl -w

use strict;
my $fasta = shift;
die "Gimme a fasta file!\n" unless $fasta && -s $fasta;
my $outfile = $fasta;
$outfile=~s/\.\w+$//;

my $orig_sep = $/;
$/ = ">";

open (IN,$fasta);
open (OUT,">$outfile.fastq");

while (my $record=<IN>){
	chomp($record);
	next unless $record;
	my @lines = split("\n",$record);
	my $id = shift(@lines);
	my $seq=join('',@lines);;
	$seq=~s/\s+$//;
	next if length($seq)<100;

	my @array = split("",$seq);
	my $qual_str;
	foreach my $a (@array){
		 if ($a eq 'N'){
			$qual_str .= '#'; 
		}else{
			$qual_str .= '=';
		}
	}

	##my $qual_str = 'C' x (length($seq)-1);
	print OUT "@".$id."\n";
	print OUT uc($seq)."\n";
	print OUT "+".$id."\n";
	print OUT $qual_str."\n";
}

close IN;
close OUT;
