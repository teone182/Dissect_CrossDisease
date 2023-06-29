#!/usr/bin/perl

#Usage: GeneRiskScore_WithLinearConversion.pl input(CADD Weight file) Outfile
#Uses CADD_raw (v1.3) scores#####


use warnings;
use strict; 

open( IN,  "$ARGV[0]" );
open( OUT, ">>$ARGV[1]" );

my @data_row;
my $Variant_ID;
my $Variant_CADD;

while ( my $line = <IN> ) {
    my @data_row;
    chomp $line;

           @data_row = split( "\t", $line );        
	   $Variant_ID = $data_row[0];
	   $Variant_CADD = $data_row[1];

           my $NEW_CADD;
           my $OldMax = "35.788538"; 
           my $OldMin = "-7.535037";
           my $NewMax = "1";
           my $NewMin = "0";


                   $NEW_CADD = ((($Variant_CADD-$OldMin)*($NewMax-$NewMin)) / ($OldMax-$OldMin)) + $NewMin;
                   print OUT "$Variant_ID\t$NEW_CADD\n";                   

} 
               

close IN;
close OUT;



