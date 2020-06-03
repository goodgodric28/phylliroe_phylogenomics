#!/usr/bin/perl

%TF_length = ();

open FASTA, "$ARGV[0]";

$TF = "";
$length = 0;

while(<FASTA>) {
    chomp;
    $line = $_;
    if($line =~ />(.*)/) {
	if($TF ne "") {
	    $TF_length{$TF} = $length;
	    $length = 0;
	}
	$TF = $1;
    } else {
	$length += length($line);
    }
}
$TF_length{$TF} = $length;
close FASTA;

$numTFs = 0;
$totalbases = 0;

open TFLENGTHS, ">tflengths";

foreach $key (sort { $TF_length{$b} <=> $TF_length{$a} } keys %TF_length ) {
    $numTFs++;
    $length = $TF_length{$key};
    $totalbases += $length;
    print TFLENGTHS $length . "\n";
}

close TFLENGTHS;

print "number of TFs: $numTFs\n";
print "total bases: $totalbases\n";

$N50 = 0;
$L50 = 0;
foreach $key (sort { $TF_length{$b} <=> $TF_length{$a} } keys %TF_length ) {
    $length = $TF_length{$key};
    $N50 += $length;
    $L50++;
    if($N50 > $totalbases / 2) {
	print "N50: " . $length . "\n";
	print "L50: " . $L50 . "\n";
	last;
    }
}

	
    
    
    
	
	

