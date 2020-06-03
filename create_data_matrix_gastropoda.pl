#!/ur/bin/perl                                                                

use lib '/PATH/TO/hamstr.v13.2.6/lib';
use lib '/PATH/TO/hamstr.v13.2.6/lib/Bio';
use lib '/PATH/TO/hamstr.v13.2.6/lib/XML';
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;

#
# codon2aa
#
# A subroutine to translate a DNA 3-character codon to an amino acid
#   Version 3, using hash lookup

sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TCN' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    #'TAA' => '_',    # Stop
    #'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    #'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CTN' => 'L',
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CCN' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'CGN' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'ACN' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GTN' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GCN' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    'GGN' => 'G',    # Glycine
    'TGA' => 'U',    # Selenocysteine
    'TAA' => 'U'     # Selenocysteine
	);

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    } else {
	print "Bad codon \"$codon\"!!\n";
	print "protchar: $protchar\n";
	print "cdspointer: $cdspointer  length of cdschars: " . scalar(@cdschars) . "\n";
	print "i: $i  header: $header\n";
	print "cds sequence: $cdsalnseq\n";
	print "skipahead: $skipahead  numlowercase: $numlowercase  amounttoadd: $amounttoadd\n";
	exit;
    }
}

$consensus = 1;
$database = "";
$outputdir = "";

if(@ARGV != 2) {
    print "Usage: perl create_data_matrix_gastropoda.pl [con|rep] [outdir]\n";
    exit;
} else {
    if($ARGV[0] eq "con") {
	$consensus = 1;
    } elsif($ARGV[0] eq "rep") {
	$consensus = 0;
    } else {
	print "First argument must be \"con\" (consensus) or \"rep\" (representative)\n";
	exit;
    }
    $outputdir = $ARGV[1];
    if(-e $outputdir) {
	print "Specified output directory $outputdir already exists! Please remove it before continuing.\n";
	exit;
    } else {
	# create output directory
	`mkdir $outputdir`;
    }
}

@hamstr_output_files = ();

if(!-e "hamstr_output_files") {
    print "Please provide a list of HaMStR output files (*.out) you would like inclded in your data matrix in a file named hamstr_output_files. Exiting...\n";
    exit;
}

$include_aplysia_ref = 0;

open TAXA, "hamstr_output_files";
while(<TAXA>) {
    chomp;
    $taxon_file = $_;
    if($taxon_file !~ /.*\.out/ && $taxon_file !~ /Aplysia_californica_reference/) {
	print "HaMStR output file name must end in .out! exiting...\n";
	exit;
    }
    if($taxon_file =~ /Aplysia_californica_reference/) {
	$include_aplysia_ref = 1;
    } else { # check for existence of file
	unless(-e $taxon_file) {
	    print "HaMStR output file not found: $taxon_file\n";
	    print "exiting...\n";
	    exit;
	}
	$cdsfile = "";
	if($taxon_file =~ /(.*)\.out/) {
	    $cdsfile = $1 . "_cds.out";
	}
	unless(-e $cdsfile) {
	    print "HaMStR output file not found: $cdsfile\n";
	    print "exiting...\n";
	    exit;
	}	
    }
    push(@hamstr_output_files, $taxon_file);
}
close TAXA;

if(@hamstr_output_files < 1) {
    print "Need at least one HaMStR output file to continue! Please check the contents of hamstr_output_files. exiting...\n";
    exit;
}

print "Command is: create_data_matrix_gastropoda.pl ";
foreach $arg (@ARGV) {
    print $arg . " ";
}
print "\n";

%protheader = ();
%protseq = ();

print "Reading HaMStR output files... ";

foreach $filename (@hamstr_output_files) {

    open PROT, "$filename";
    
    $group = "";
    $header = "";
    $sequence = "";
    
    while(<PROT>) {
	chomp;
	$line = $_;

	# headers look like:
	# PROT: A5GZV0|Aplysia_californica|Euplana_gracilis|c34363_g1_i1|1|VNEAKKRI ...
	#  CDS; A5GZV0|Aplysia_californica|Euplana_gracilis|c34363_g1_i1|1|GTCAACGA ...

	if($line =~ /(.*?)(\|.*)\|(.*)$/) {
	    $group = $1;
	    $header = $1 . $2;
	    #print "header is: $header\n";
	    $sequence = $3;
	    
	    if(defined($protheader{$group})) {
		$protheader{$group} = $protheader{$group} . "," . $header;
		$protseq{$group} = $protseq{$group} . "," . $sequence;
	    } else {
		$protheader{$group} = $header;
		$protseq{$group} = $sequence;
	    }
	} else {
	    print "Erroneous line in $filename: $line\n";
	    #exit;
	}
    }
    close PROT;
}
	
print "done.\n";
print "Organizing sequences into orthologous groups... ";

@gastropoda_ogs = ();
open OGS, "gastropoda_ogs";
while(<OGS>) {
    chomp;
    $gastro_og = $_;
    push(@gastropoda_ogs, $gastro_og);
}
close OGS;

foreach $i (@gastropoda_ogs) { # loop through all genes
    $ogfilename = "$outputdir/" . $i . ".fa";
    open OG, ">$ogfilename";

    if($include_aplysia_ref == 1) {
	$aplysia_file = "";
	$aplysia_file = "/PATH/TO/hamstr.v13.2.6/gastropoda_50_aplysia_californica/gastropoda_50_aplysia_californica/fa_dir/$i.fa2";

	open APLYSIA, "$aplysia_file";
	
	$header = "";
	$sequence = "";
	%aplysia_sequences = ();
	while(<APLYSIA>) {
	    chomp;
	    $line = $_;
	    
	    if($line =~ />.*Aplysia_californica.*/) {
		$header = $line;
		$sequence = <APLYSIA>;
		chomp($sequence);
		$aplysia_sequences{$header} = $sequence;
	    }
	}
	if($consensus == 1) { # simply print all sequences
	    foreach $aplysia_header (sort keys %aplysia_sequences) {
		$aplysia_sequence = $aplysia_sequences{$aplysia_header};
		print OG $aplysia_header . "\n" . $aplysia_sequence . "\n";
	    }
	} else { # print the longest sequence
	    $longest_sequence_header = "";
	    $longest_sequence_length = 0;
	    foreach $aplysia_header (sort keys %aplysia_sequences) {
		$aplysia_sequence = $aplysia_sequences{$aplysia_header};
		$aplysia_sequence_length = length($aplysia_sequence);
		if($aplysia_sequence_length > $longest_sequence_length) {
		    $longest_sequence_header = $aplysia_header;
		    $longest_sequence_length = $aplysia_sequence_length;
		}
	    }
	    print OG $longest_sequence_header . "\n" . $aplysia_sequences{$longest_sequence_header} . "\n";
	}
	close APLYSIA;
    }

    $group = $i;

    if(defined($protheader{$group})) {
	@headers = split(/,/, $protheader{$group});
	@sequences = split(/,/, $protseq{$group});
	for($j = 0; $j < @headers; $j++) {
	    $header = $headers[$j];
	    $sequence = $sequences[$j];
	    print OG ">" . $header . "\n" . $sequence . "\n";
	}
    }
    close OG;
}
print "done.\n";

# create CDS directory
`mkdir $outputdir/CDS`;

print "Reading HaMStR CDS output files... ";

$cdsheader = ();
$cdsseq = ();

foreach $filename (@hamstr_output_files) {
    if($filename =~ /(.*)\.out/) {
	$filename = $1 . "_cds.out";
    }

    open CDS, "$filename";
    
    $group = "";
    $header = "";
    $sequence = "";
    
    while(<CDS>) {
	chomp;
	$line = $_;
	
	if($line =~ /(.*?)(\|.*)\|(.*)$/) {
	    $group = $1;
	    $header = $1 . $2;
	    $sequence = $3;
	    
	    if(defined($cdsheader{$group})) {
		$cdsheader{$group} = $cdsheader{$group} . "," . $header;
		$cdsseq{$group} = $cdsseq{$group} . "," . $sequence;
	    } else {
		$cdsheader{$group} = $header;
		$cdsseq{$group} = $sequence;
	    }
	} else {
	    print "Erroneous line in $filename: $line\n";
	    #exit;
	}
    }
    close CDS;
}

print "done.\n";
print "Organizing CDS sequences into orthologous groups... ";

foreach $i (@gastropoda_ogs) { # loop through all genes
    $ogfilename = "$outputdir/CDS/" . $i . ".fa";
    open OG, ">$ogfilename";

    # we are currently not writing out Aplysia CDS sequences
    #if($include_aplysia_ref == 1) {
    # ...      
    #}

    $group = $i;
    
    if(defined($cdsheader{$group})) {
	@headers = split(/,/, $cdsheader{$group});
	@sequences = split(/,/, $cdsseq{$group});
	for($j = 0; $j < @headers; $j++) {
	    $header = $headers[$j];
	    $sequence = $sequences[$j];
	    print OG ">" . $header . "\n" . $sequence . "\n";
	}
    }
    close OG;
}

print "done.\n";

print "Replacing '*' characters with 'X'...";
`perl -pi -e 's/\\*/X/g' $outputdir/*.fa`;
print "done.\n";

print "Aligning orthologous group sequences with MAFFT, this may take a few minutes... ";

# create alignment directory
`mkdir $outputdir/mafft_aligned`;

# copy in saved alignments
#`cp -r saved_alignments/* $outputdir/mafft_aligned`;

#=pod

foreach $i (@gastropoda_ogs) { # loop through all genes
    
    $in = $outputdir . "/" . $i . ".fa";
    $out = $outputdir . "/mafft_aligned/$i.aln"; # final output file
    
    $numseqs = `grep -c ">" $in`;
    chomp($numseqs);
    
    if($numseqs eq "1") {
        `cp $in $out`;
    } else {
	open INPUT, "$in";
	@headers = ();
	@sequences = ();
	$sequence = "";
	while(<INPUT>) {
	    chomp;
	    $line = $_;
	    if($line =~ />.*/) {
		if($sequence ne "") {
		    push(@sequences, $sequence);
		    $sequence = "";
		}
		push(@headers, $line);
	    } else {
		$sequence .= $line;
	    }
	}
	if($sequence ne "") {
	    push(@sequences, $sequence);
	}
	close INPUT;

	# create temp dir in which to work
	$workingdir = $outputdir . "/mafft_aligned/$i";
	`mkdir $workingdir`;

	$refheader = $headers[0];
	$refsequence = $sequences[0];

	# create aplysia file
	$aplysiafile = $workingdir . "/aplysia";
	open APLYSIA, ">$aplysiafile";
	print APLYSIA $refheader . "\n";
	print APLYSIA $refsequence . "\n";
	close APLYSIA;

	# create TFs file
	$tfsfile = $workingdir . "/TFs";
	open TFS, ">$tfsfile";

	for($j = 1; $j < @headers; $j++) {
	    $header = $headers[$j];
	    $sequence = $sequences[$j];
	    print TFS $header . "\n";
	    print TFS $sequence . "\n";
	}
	close TFS;

	# align all fragments with Aplysia
	#`mafft --thread -1 --addfragments $tfsfile $bombyxfile > $out 2>> $outputdir/mafft_aligned/mafft.log`;
	$mafftoutput = `mafft --aamatrix /PATH/TO/hamstr.v13.2.6/custom_blast_matrices/GASTRO50 --anysymbol --addfragments $tfsfile $aplysiafile > $out 2> $outputdir/mafft_aligned/$i.mafft.log`;

	# if alignment file is zero-size, exit with error
	if((-s "$out") == 0) {
	    print "alignment output for OG $i is empty! removing output file...\n";
	    `rm $out`;
	    #print "alignment output for OG $i is empty! exiting...\n";
	    #exit;
	}
    }
}
#=cut

print "done.\n";
print "Converting amino acid alignments to nucleotide alignments... ";

$ogs_skipped = 0;

foreach $i (@gastropoda_ogs) { # loop through all genes

    if(-e "$outputdir/mafft_aligned/$i.aln") {
	%protseq = ();
	%cdsseq = ();
	open PROTALN, "$outputdir/mafft_aligned/$i.aln";
	open CDSSEQ, "$outputdir/CDS/$i.fa";
	$cdsalnfile = $outputdir . "/mafft_aligned/" . $i . "_cds.aln";
	open CDSALN, ">$cdsalnfile";

	@headerlist = ();
	
	$header = "";
	$sequence = "";
	
	while(<PROTALN>) {
	    chomp;
	    $line = $_;
	    
	    if($line =~ />.*/) {
		if($header ne "") {
		    $protseq{$header} = $sequence;
		    push(@headerlist, $header);
		}
		$header = $line;
		$sequence = "";
	    } else {
		$sequence .= $line;
	    }
	}
	if($header ne "") {
	    $protseq{$header} = $sequence;
	    push(@headerlist, $header);
	}
	close PROTALN;
	
	$header = "";
	$sequence = "";
	
	while(<CDSSEQ>) {
	    chomp;
	    $line = $_;
	    
	    if($line =~ />.*/) {
		if($header ne "") {
		    $cdsseq{$header} = $sequence;
		}
		$header = $line;
		$sequence = "";
	    } else {
		$sequence .= $line;
	    }
	}
	if($header ne "") {
	    $cdsseq{$header} = $sequence;
	}
	close CDSSEQ;
	
	for($h = 0; $h < @headerlist; $h++) {
	    $aplysia = 0;
	    $header = $headerlist[$h];
	    
	    if($header =~ /Aplysia_californica/) { # this is and Aplysia sequence JAG - No, it's not. Need to revisit this.
		$aplysia = 1;
	    }
	    
	    $cdsalnseq = "";
	    $protsequence = $protseq{$header};
	    $cdssequence = "";
	    
	    if(defined($cdsseq{$header}) && $cdsseq{$header} ne "") {
		print CDSALN $header . "\n";
		$cdssequence = $cdsseq{$header};
		@protchars = split(undef, $protsequence);
		@cdschars = split(undef, $cdssequence);
		
		$prot_seq_without_gaps = $protsequence;
		$prot_seq_without_gaps =~ s/\-//g;
		$readingframe = 0;
		
		if(length($cdssequence) < (length($prot_seq_without_gaps) * 3)) {
		    print "not enough CDS characters! header: $header\n";
		    print "length of cds sequence: " . length($cdssequence) . "  length of prot sequence: " . length($prot_seq_without_gaps) . "\nexiting...\n";
		    exit;
		}
		
		$cdspointer = 0;
		for($j = 0; $j < @protchars; $j++) {
		    $protchar = $protchars[$j];
		    if($protchar eq "-") {
			$cdsalnseq .= "---";
		    } elsif($protchar eq "X" || $protchar eq "x") {
			print "protchar is $protchar\n";
			# if this is the last amino acid, simply grab the corresponding codon
			#print "j: $j\n";
			#print "length of protsequence: " . length($protsequence) . "\n";
			#print "length of protchars: " . @protchars . "\n";
			#print "protsequence: " . $protsequence . "\n";
			$allgaps = 1;
			for($k = $j + 1; $k < @protchars; $k++) {
			    if($protchars[$k] ne "-") {
				$allgaps = 0;
				last;
			    }
			}
			
			if($allgaps == 1) { # the x character was the last amino acid
			    $cdsalnseq .= $cdschars[$cdspointer] . $cdschars[$cdspointer+1] . $cdschars[$cdspointer+2];
			    $cdspointer += 3;
			} else {
			    $tempcodon = $cdschars[$cdspointer] . $cdschars[$cdspointer+1] . $cdschars[$cdspointer+2];
			    print "tempcodon is: $tempcodon\n";
			    
			    # we need to determine if this represents a stop codon
			    if($tempcodon eq "TAG" || $tempcodon eq "TGA" || $tempcodon eq "TAA") {
				$cdsalnseq .= "NNN";
				$cdspointer += 3;
				
				# check for the case where the stop codon is followed by lowercase letters, and skip over those
				# (as in 413055 TGAt--)
				
				$numlowercase = 0;
				$temppointer = $cdspointer;
				while($cdschars[$temppointer] =~ /[acgtnryswkmbdhv]/) { # allow nucleotide ambiguity codes as well
				    $numlowercase++;
				    $temppointer++;
				}
				
				if($numlowercase > 0) {
				    print "found a case where stop codon is followed by lowercase\n";
				    print "i: $i  cdsseq: $cdssequence\n";
				}
				
				# numlowercase needs to be rounded up to a multiple of 3 (codon length)
				$remainder = $numlowercase % 3;
				$amounttoadd = 0;
				if($remainder > 0) { # we need to pad out this codon
				    $amounttoadd = 3 - $remainder;
				}
				
				$skipahead = $numlowercase + $amounttoadd;
				
				$cdspointer += $skipahead;
				
			    } else {
				
				$cdsalnseq .= "NNN";
				# we need to figure out how many codons to skip ahead
				$numlowercase = 0;
				
				# we expect at least one lowercase character introduced by HaMStR
				if($cdschars[$cdspointer] =~ /[acgtnryswkmbdhv]/) { # allow nucleotide ambiguity codes as well
				    $numlowercase++;
				} else {
				    print "expected a lowercase nucleotide; instead, got " . $cdschars[$cdspointer] . "\n";
				    print "i is: $i\n";
				    print "header is: $header\n";
				    print "cds sequence is: $cdssequence\n";
				    exit;
				}
				
				# get the rest of the lowercase characters
				$temppointer = $cdspointer+1;
				while($cdschars[$temppointer] =~ /[acgtnryswkmbdhv]/) { # allow nucleotide ambiguity codes as well
				    $numlowercase++;
				    $temppointer++;
				}
				
				# numlowercase needs to be rounded up to a multiple of 3 (codon length)
				$remainder = $numlowercase % 3;
				$amounttoadd = 0;
				if($remainder > 0) { # we need to pad out this codon
				    $amounttoadd = 3 - $remainder;
				}
				
				$skipahead = $numlowercase + $amounttoadd;
				
				$cdspointer += $skipahead;
			    }
			}
		    } elsif($protchar =~ /[A-Z]/) {
			$codon = $cdschars[$cdspointer] . $cdschars[$cdspointer+1] . $cdschars[$cdspointer+2];
			
			# make sure this codon is valid for the amino acid
			$aa = codon2aa($codon);
			
			if($aa ne $protchar) {
			    print "codon $codon does not code for $protchar! it codes for $aa...\n";
			    print "cdspointer: $cdspointer  length of cdschars: " . scalar(@cdschars) . "\n";
			    print "i: $i  header: $header\n";
			    print "cds sequence: $cdsalnseq\n";
			    print "skipahead: $skipahead  numlowercase: $numlowercase  amounttoadd: $amounttoadd\n";
			    exit;
			} else {
			    $cdsalnseq .= $codon;
			    $cdspointer += 3;
			}
		    } else {
			print "unknown protein character: $protchar\n";
			exit;
		    }
		}
		
		# cds pointer should now be pointing at @cdschars; if not, this may be a problem
		if($cdspointer != scalar(@cdschars)) {
		    print "cdspointer not equal to length of cdschars!\n";
		    print "cdspointer: $cdspointer  length of cdschars: " . scalar(@cdschars) . "\n";
		    print "i: $i  header: $header\n";
		    print "cds sequence: $cdsalnseq\n";
		    print "skipahead: $skipahead  numlowercase: $numlowercase  amounttoadd: $amounttoadd\n";
		    exit;
		}
		
		@sixties = unpack("A60" x (length($cdsalnseq)/60), $cdsalnseq);
		for $cdsalnline(@sixties) {
		    print CDSALN $cdsalnline . "\n";
		}
		# print the last line, if there is one
		if(length($cdsalnseq) % 60 != 0) {
		    print CDSALN substr($cdsalnseq, -(length($cdsalnseq) % 60))."\n";
		}
	    } else {
		print "CDS sequence not found for header: $header\n";
		#exit; # continue on, assuming these few missing sequences are due to HaMStR bugs and/or concurrent filesystem write mistakes
                       # or because we don't have CDS for Aplysia!
	    }
	}
	close CDSALN;
    } else {
	$ogs_skipped++;
    }
}
print "done.\n";

if($consensus == 1) {

    print "Preparing directories for consensus method... ";

    foreach $i (@gastropoda_ogs) { # loop through all genes
	%cdsalnseq = ();
	$cdsalnfile = $outputdir . "/mafft_aligned/" . $i . "_cds.aln";
	if(-e $cdsalnfile) {
	    open CDSALN, "$cdsalnfile";
	    
	    $header = "";
	    $sequence = "";
	    
	    while(<CDSALN>) {
		chomp;
		$line = $_;
		
		if($line =~ />.*/) {
		    if($header ne "") {
			$cdsalnseq{$header} = $sequence . "\n";
		    }
		    $header = $line;
		    $sequence = "";
		} else {
		    $sequence .= $line;
		}
	    }
	    if($header ne "") {
		$cdsalnseq{$header} = $sequence . "\n";
	    }
	    close CDSALN;
	    
	    $dir = $outputdir . "/mafft_aligned/consensus/$i";
	    `mkdir -p $dir`;
	    
	    foreach $header (keys %cdsalnseq) {
		$taxon = "";
		
		if($header =~ /.*?\|.*?\|(.*?)\|.*/) {
		    $taxon = $1;
		} else { # is a reference sequence of some kind
		    print "taxon unknown! exiting...\n";
		    exit;
		}
		$sequence = $cdsalnseq{$header};
		
		$file = $dir . "/" . $i . "_$taxon\_cds.aln";	
		open FILE, ">>$file";
		print FILE $header . "\n" . $sequence;
		close FILE;
	    }
	}
    }
    
    print "done.\n";
    print "Creating consensus sequences... ";
    
    $length = 0;
    $ambi = 0;

    foreach $i (@gastropoda_ogs) { # loop through all genes

	$cdsalnfile = $outputdir . "/mafft_aligned/" . $i . "_cds.aln";
	if(-e $cdsalnfile) {
	    
	    $path = $outputdir . "/mafft_aligned/consensus/$i";
	    @alnfiles = <$path/*.aln>; 
	    
	    foreach $alnfile (@alnfiles) {
		
		$taxon = "";
		if($alnfile =~ /.*_(.*)_cds.aln/) {
		    $taxon = $1;
		} else {
		    print "Could not extract taxon from alnfile $alnfile - exiting...\n";
		    exit;
		}
		
		$out = $path . "/" . $i . "_con_cds.aln";
		open OUT, ">>$out";
		
		$alignio = new Bio::AlignIO(-format => 'fasta', -file => "$alnfile");
		$aln = $alignio->next_aln();
		if(defined($aln)) {
		    for $seq ($aln->each_seq) {
			$seq->alphabet('dna');
		    }
		    $con = uc($aln->consensus_iupac());
		    
		    @conchars = split(undef, $con);
		    foreach $char(@conchars) {
			if($char !~ /[ACGT\-]/) {
			    $ambi++;
			    $length++;
			} elsif($char !~ /\-/) {
			    $length++;
			}
		    }
		    
		    print OUT ">$taxon\_CON\n";
		    @sixties = unpack("A60" x (length($con)/60), $con);
		    for $line(@sixties) {
			print OUT $line . "\n";
		    }
		    # print the last line, if there is one
		    if(length($con) % 60 != 0) {
			print OUT substr($con, -(length($con) % 60)) . "\n";
		    }
		}
		
		close OUT;
	    }
	}
    }
    print "done.\n";

    print "number of characters: $length\n";
    print "number of ambiguous characters: $ambi\n";
    $ambifrac = $ambi / $length;
    print "fraction of ambiguous characters: $ambifrac\n";
}

exit;

print "Concatenating genes and adding gaps for missing data... ";

# create concatenated directory
`mkdir $outputdir/mafft_aligned/concatenated`;

$first_file = 1;
%sequence_ids = ();
%sequences = ();
%sequence_lengths = ();
%master_sequences = ();

# open final combined file
#$number_of_genes = (($last_og - $first_og) + 1) - $ogs_skipped;
$finaldatamatrix = $outputdir . "/mafft_aligned/concatenated/" . $outputdir . "_" . $number_of_genes . "_genes.fa";
print "\n\n outputdir: $outputdir\n\n finaldatamatrix: $finaldatamatrix\n\n";
open FINALDATAMATRIX, ">$finaldatamatrix";

foreach $i (@gastropoda_ogs) { # loop through all genes

    if($consensus == 1) {
	$file = $outputdir . "/mafft_aligned/consensus/" . $i . "/" . $i . "_con_cds.aln";
    } else { # must be using -representative
	$file = $outputdir . "/mafft_aligned/" . $i . "_cds.aln";
    }

    if(-e $file) { # some files may be missing due to OGs having been removed from database

	$FASTA_in = Bio::SeqIO->new(-file => $file);
	while($fasta_seq = $FASTA_in->next_seq()) {
	    
	    $taxon = "";
	    if($consensus == 1) {
		if($fasta_seq->primary_id =~ /(.*?)_CON/) {
		    $taxon = $1;
		} else {
		    print "Could not extract taxon from consensus sequence identifier! exiting...\n";
		    exit;
		}
	    } else {
		if($fasta_seq->primary_id =~ /.*?\|.*?\|(.*?)\|.*/) {
		    $taxon = $1;
		} else { # is a reference sequence of some kind
		    if($header =~ /BGIBMGA/) { 
			$taxon = "Bombyx_mori";
		    } elsif($header =~ /KGM_/) {
			$taxon = "Danaus_plexippus";
		    } elsif($header =~ /HMEL/) {
			$taxon = "Heliconius_melpomene";
		    } elsif($header =~ /Unigene/ || $header =~ /Px/) {
			$taxon = "Plutella_xylostella";
		    } elsif($header =~ /Msex/) {
			$taxon = "Manduca_sexta";
		    } else {
			print "taxon unknown! exiting...\n";
			exit;
		    }
		}
	    }
	    
	    # push the sequence id
	    push( @{$sequence_ids{$taxon}}, $file );
	    
	    # push the sequence data
	    push( @{$sequences{$taxon}}, $fasta_seq->seq );
	    
	    # push the sequence length
	    $sequence_lengths{$file} = length($fasta_seq->seq);
	}
    }
}

$genelengthsfile = $outputdir . "/mafft_aligned/concatenated/gene_lengths";
open GENELENGTHS, ">$genelengthsfile";
$didfirst = 0;

# loop through the hash of sequence ids
foreach $seq_id (sort keys %sequence_ids) {

    foreach $i (@gastropoda_ogs) { # loop through all genes
	
	if($consensus == 1) {
	    $file = $outputdir . "/mafft_aligned/consensus/" . $i . "/" . $i . "_con_cds.aln";
	} else { # must be using -representative
	    $file = $outputdir . "/mafft_aligned/" . $i . "_cds.aln";
	}

	if(-e $file) {
	    
	    if($didfirst == 0) {
		$genelength = $sequence_lengths{$file};
		print GENELENGTHS $genelength . "\n";
		$didfirst = 1;
	    }
	    
	    $file_found = 0;
	    $counter = 0;
	    
	    foreach $seq_file (@{$sequence_ids{$seq_id}}) {
		if($seq_file eq $file) {
		    $file_found = 1;
		    
		    # append the sequence data
		    $master_sequences{$seq_id} = $master_sequences{$seq_id} . ${$sequences{$seq_id}}[$counter];
		    last;
		}
		$counter++;
	    }
	    
	    if($file_found == 0) { # file was not found, so append some dashes
		my $dashes = "";
		for(my $x = 0; $x < $sequence_lengths{$file}; $x++) {
		    $dashes = $dashes."-";
		}
		$master_sequences{$seq_id} = $master_sequences{$seq_id} . $dashes;
	    }
	}
    }
}
close GENELENGTHS;

# loop through the hash of sequence ids
foreach $seq_id (sort keys %sequence_ids) {

    # retrieve sequence information
    $combined_sequence = $master_sequences{$seq_id};

    $seq_id =~ s/,//g; # remove any commas

    print FINALDATAMATRIX ">" . $seq_id . "\n";

    # now, output the sequences in lines of sixty characters
    @sixties = unpack("A60" x (length($combined_sequence)/60), $combined_sequence);

    foreach $line (@sixties) {
        print FINALDATAMATRIX $line . "\n";
    }

    # print the last line, if there is one
    if(length($combined_sequence) % 60 != 0) {
        print FINALDATAMATRIX substr($combined_sequence,
                            -(length($combined_sequence) % 60)) . "\n";

    }
}

close FINALDATAMATRIX;

print "done.\n";

print "Length of alignment before removing columns is... ";

$alignmentlength = 0;

$alignio = new Bio::AlignIO(-format => 'fasta', -file => "$finaldatamatrix");
$aln = $alignio->next_aln();
$printed_alignment_length = 0;
if(defined($aln)) {
    for $seq ($aln->each_seq) {
        $seq->alphabet('dna');
	if($printed_alignment_length == 0) {
	    $alignmentlength = length($seq->seq);
	    print "$alignmentlength nucleotides.\n";
	    $printed_alignment_length = 1;
	}
    }
}

print "Removing phylogenetically uninformative positions... ";

%less_than_four_non_gap_positions = ();

@master_char_array = ();
foreach $seqid (sort keys %sequence_ids) {
    $sequence = $master_sequences{$seqid};
    @chars = split(//, $sequence);
    push @master_char_array, [ @chars ];
}

for($i = 0; $i < $alignmentlength; $i++) {
    $gapcount = 0;
    $nongapcount = 0;
    for($j = 0; $j < @hamstr_output_files; $j++) {
	if($master_char_array[$j][$i] eq "-") {
	    $gapcount++;
	} else {
	    $nongapcount++;
	}
    }
    if($nongapcount < 4) { # this position has less than four non-gap chars
	$less_than_four_non_gap_positions{$i} = 1;
    }
}

# open final combined file with phylogenetically uninformative positions removed
$finaldatamatrix2 = $outputdir . "/mafft_aligned/concatenated/" . $outputdir . "_" . $number_of_genes . "_genes_cols_removed.fa";
open FINALDATAMATRIX, ">$finaldatamatrix2";
foreach $seqid (sort keys %sequence_ids) {
    $sequence = $master_sequences{$seqid};
    @chars = split(//, $sequence);
    print FINALDATAMATRIX ">$seqid\n";
    for ($i = 0; $i < @chars; $i++) {
	if(!defined($less_than_four_non_gap_positions{$i})) { # print this character
	    print FINALDATAMATRIX $chars[$i];
	}
    }
    print FINALDATAMATRIX "\n";
}
close FINALDATAMATRIX;

print "done.\n";

print "Converting alignment to NEXUS format... ";
$command = "clustalw2 -infile=" . $finaldatamatrix2 . " -convert -output=nexus -outfile=" . $outputdir . "/mafft_aligned/concatenated/" . $outputdir . "_" . $number_of_genes . "_genes_cols_removed.nex";
$output = `$command`;
print "done.\n";

print "Length of alignment after removing columns is... ";

$alignmentlength = 0;
%nongaptotals = ();
$taxontotalnongap = 0;
$totalnongap = 0;

$alignio = new Bio::AlignIO(-format => 'fasta', -file => "$finaldatamatrix2");
$aln = $alignio->next_aln();
$printed_alignment_length = 0;
if(defined($aln)) {
    for $seq ($aln->each_seq) {
        $seq->alphabet('dna');
        @chars = split(//,$seq->seq);
        foreach $char (@chars) {
            if($char ne '-') {
                $taxontotalnongap++;
            }
        }
	if($printed_alignment_length == 0) {
	    $alignmentlength = length($seq->seq);
	    print "$alignmentlength nucleotides.\n";
	    $printed_alignment_length = 1;
	}

        $id = $seq->id;
	$nongaptotals{$id} = $taxontotalnongap;
	$totalnongap += $taxontotalnongap;

        $taxontotalnongap = 0;
    }
}

$total_possible_genes = @hamstr_output_files * $number_of_genes;

=pod
print "Calculating matrix completeness by gene... ";
$genes_present = 0;
if($consensus == 1) {
    $condir = $outputdir . "/mafft_aligned/consensus";
    $command = "find $condir -name '*con*' -exec grep -c '>' {} \\; | awk '{s+=\$0} END {print s}'";
    print "command is: $command\n";
    $genecount = `$command`;
    if($genecount =~ /([0-9]*)/) {
	$genecount = $1;
    }
    $genes_present = $genecount;

    # add bombyx
    $genes_present += $number_of_genes;
} else { # -representative
    $alndir = $outputdir . "/mafft_aligned";
    $command = "find $alndir -maxdepth 1 -name '*_cds.aln' -exec grep -c '>' {} \\; | awk '{s+=\$0} END {print s}'";
    $genecount = `$command`;
    if($genecount =~ /([0-9]*)/) {
	$genecount = $1;
    }
    $genes_present = $genecount;
    if($include_bombyx_ref == 1) {
	# add bombyx
	$genes_present += $number_of_genes;
    }
}
print "done.\n";
=cut

print "Calculating matrix completeness by position... ";
$total_possible_positions = @hamstr_output_files * $alignmentlength;
$completeness_by_position = sprintf("%.2f", ($totalnongap / $total_possible_positions) * 100);
print "done.\n";

# output completeness statistics
print "Outputting completeness statistics... ";
open COMPLETENESS, ">$outputdir/mafft_aligned/concatenated/completeness_stats";
#print COMPLETENESS $genes_present . " genes present out of " . $total_possible_genes . " total possible genes = " . sprintf("%.2f", ($genes_present / $total_possible_genes) * 100) . "\% completeness by gene\n\n";
print COMPLETENESS $totalnongap . " positions present out of " . $total_possible_positions . " total possible positions = " . $completeness_by_position . "\% completeness by position\n\n";

foreach $taxon (sort keys %nongaptotals) {
    print COMPLETENESS "  taxon: $taxon  non-gap characters: " . $nongaptotals{$taxon} . "\n";
}

close COMPLETENESS;

print "done.\n";
print "Creating gene charsets... ";

@genelengths = ();
open GENELENGTHS, "$genelengthsfile";
while(<GENELENGTHS>) {
    chomp;
    $genelength = $_;
    push(@genelengths, $genelength);
}
close GENELENGTHS;

$genecharsetfile = $outputdir . "/mafft_aligned/concatenated/" . $outputdir . "_" . $number_of_genes . "_genes--gene_charsets";
open CHARSET, ">$genecharsetfile";

$genecounter = 1;
$positioncounter = 1;
foreach $genelength (@genelengths) {
    $endextent = ($positioncounter+$genelength) - 1;
    $charsetstring = "charset gene_" . $genecounter . " = " . $positioncounter . "-" . $endextent . ";";
    print CHARSET $charsetstring . "\n";
    $positioncounter += $genelength;
    $genecounter++;
}

close CHARSET;

print "done.\n";

print "Output files can be found in " . $outputdir . "/mafft_aligned/concatenated\n";
