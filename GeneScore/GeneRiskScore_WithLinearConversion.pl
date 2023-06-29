#!/usr/bin/perl

#Usage: GeneRiskScore_WithLinearConversion.pl input(.meta file) Outfile
#Uses CADD_raw (v1.3) scores#####


use warnings;
use strict; 

open( IN,  "$ARGV[0]" );
open( OUT, ">>$ARGV[1]" );

$" = "\t";

while ( my $line = <IN> ) {
    my @genotypes_row;
    chomp $line;
    if ( $line =~ m/^Chr/ ) {
        print OUT "$line\n";
    }
    else {
           my $Variant_CADD;

           @genotypes_row = split( "\t", $line );        
	   my $Variant_Chrom = $genotypes_row[0];
	   my $Variant_Pos = $genotypes_row[1];
	   my $Variant_Gene = $genotypes_row[6];

           if ($genotypes_row[10] eq ".") {
               $Variant_CADD = "-7.535037";
            }
                else {
                $Variant_CADD = $genotypes_row[10]        
                }

           my $new_elem;
           my $new_elem_CADD;
           my $OldMax = "35.788538"; 
           my $OldMin = "-7.535037";
           my $NewMax = "1";
           my $NewMin = "0";
           my $CADD_conv;

           for ( my $i = 0 ; $i <= $#genotypes_row ; $i++ ) {
           my $elem = $genotypes_row[$i];           

               #if (($elem =~ m/^[aA-zZ]/) || ($elem =~ m/^\.$/) || ($elem =~ m/^\d+\.?\d+$/) || ($elem =~ m/^\-\d+\.?\d+$/)) {
               if (($elem =~ m/^[aA-zZ]/) || ($elem =~ m/^\.$/) || ($elem =~ m/^-?[0-9]\d*(\.\d+)?$/)) {
               $new_elem_CADD = $elem;
               }
               elsif (($elem =~ m/^0\/0:/) || ($elem =~ m/^0:/)) {
               $new_elem = 0 * 1.0; 
               $new_elem_CADD = $new_elem;
               }
                
               elsif (($elem =~ m/^0\/1:/) || ($elem =~ m/^1\/0:/) || ($elem =~ m/^1:/)) {
                   $CADD_conv = ((($Variant_CADD-$OldMin)*($NewMax-$NewMin)) / ($OldMax-$OldMin)) + $NewMin;
                   $new_elem_CADD = $CADD_conv*1;
                   } 
               
               elsif ( $elem =~ m/^1\/1:/) {
                   $CADD_conv = ((($Variant_CADD-$OldMin)*($NewMax-$NewMin)) / ($OldMax-$OldMin)) + $NewMin;
                   $new_elem_CADD = $CADD_conv*2;
                   }

               elsif (($elem =~ m/^\.\/\.:/) || ($elem =~ m/^\.:/)) {
               $new_elem = 0 * 1.0;      ######0 instead of NA####
               $new_elem_CADD = $new_elem;
               }
               else { ; }

            $genotypes_row[$i] = $new_elem_CADD;
            }
            print OUT "@genotypes_row\n";
    }
}

close IN;
close OUT;




