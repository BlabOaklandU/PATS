#!/usr/bin/env perl
use strict;
use warnings "all";

if (!defined($ARGV[0]) || $ARGV[0] eq "-h" || $ARGV[0] =~ /^-?-help/) {
	print STDERR "proteinortho_singletons.pl FASTA1 FASTA2 ... FASTAN <PROTEINORTHO_OUTFILE >SINGLETON_GENES\n";
	print STDERR "Reads Proteinortho outfile and its source fasta files to determin entries which occure once only\n\n";
	exit;
}

sub convertUniprotAndNCBI {
  my $long_id = shift;
  $long_id =~ s/\|$//g;
  my @tmp = split(/\|/,$long_id); 
  if(scalar @tmp > 3){
    return pop(@tmp); # take the last column for NCBI format like patterns (e.g. gi|158333234|ref|YP_001514406.1|)
  }elsif(scalar @tmp == 3){
    return $tmp[1]; # uniprot tr|A0A0A0MPE6|A0A0A0MPE6_DANRE -> A0A0A0MPE6
  }else{
    return $long_id; # is neither ncbi nor uniprot 
  } 
}

my %present; # Genes present in the matrix

### Parse matrix, store present genes and species order
my %order;
while (<STDIN>) {
	if ($_ =~ /#\sSpecies/) {
		my @species = split(/\s+/);
		shift @species;shift @species;shift @species;shift @species;
		for (my $i = 0; $i < scalar(@species); $i++) {
			$order{$species[$i]} = $i;
		}
		next;
	}
	if ($_ =~ /^#/ || length($_) < 4) {next;}
	chomp;
	my @row = split(/\s|,/,$_);
	for (my $i = 3; $i < scalar(@row); $i++) {
		$present{$row[$i]} = 1;
	}
}

### For each fasta file, parse gene-IDs
foreach my $file (@ARGV) {
	# print "\# $file\n";
	my $pos;
	unless (defined($order{$file})) {
		my $short_version = $file;
		$short_version =~ s/^.*\///;
		unless (defined($order{$short_version})) {
			die("Species $file is not in the matrix\n");
		}
		$pos = $order{$short_version};		
	}
	else {$pos = $order{$file};}
	open(FILE,"<$file") || die("Error, could not open file $file: $!");
	while(<FILE>) {
		if ($_ =~ /^>([^\s]+)/) {
			my $id = $1;
			if (exists($present{$id})) {next;}
			if (exists($present{&convertUniprotAndNCBI($id)})) {next;} # mmseqs converts the ids automatically -> fix for -keep -p=mmseqsp
			
			# add to matrix
			print "1\t1\t0";
			for (my $i = 0; $i < scalar(keys %order); $i++) {
				print "\t";
				if ($i == $pos) {
					print $id;
				}
				else {
					print "*";
				}
			}
			print "\n";
		}
	}
	close(FILE);
}

