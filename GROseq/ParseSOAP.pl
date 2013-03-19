#!/usr/bin/perl
# Program grabs the columns from a SOAP2 run and writes them to 
# standard output ...
my $readSize = 44;
my $name = "n";

while(<STDIN>) {
	@SPL = split("\t");
	$strand     = @SPL[6];
	$chrom      = @SPL[7];
	$chromStart = @SPL[8];
	$chromEnd   = $chromStart+$readSize;
	$score      = @SPL[9];

	if($strand eq "+") {
		$FivePrimeStart = $chromStart
	}
	else {
		if($strand eq "-") {
			$FivePrimeStart = $chromEnd
		}
	}

	print "$chrom\t$FivePrimeStart\t".($FivePrimeStart+1)."\t$name\t$score\t$strand\n";
}
