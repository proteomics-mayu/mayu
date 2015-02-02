#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use FileHandle;

##################################################################
# 
# use: reverse_fasta.pl in.fasta
# out: a target decoy database with changed headers and reversed 
#      sequences
# 
##################################################################

#-------------------------------------------------------------
# PARAMETERS AND INPUT
#-------------------------------------------------------------
my $v = 1;
my $s = 1;

# added as a prefix to the ids of the proteins
my $prefix = 'rev_';

# name of a fasta protein database
my $indb = shift;

unless (defined($indb)) {
	print "  use: $0 in.fasta\n";
	exit;
}
unless ( -f $indb ) {
	print "  $indb not a file!\n";
	exit;
}

# open a output filehandle
my $outdb = 'target_decoy_' . basename($indb);

print "\n";
print "  $0\n";
print "  input fasta database: $indb\n";
print "  decoy protein id prefix: $prefix\n";
print "  output fasta target decoy database: $outdb\n";

#-------------------------------------------------------------
# START THE PROGRAM
#-------------------------------------------------------------
my $outfh = FileHandle->new();
$outfh->open( ">$outdb" ) or die $!;

# print the target protein entries
open(I, $indb) or die $!;
while(<I>) {
	my $line = $_;
	print $outfh $line;
}
close(I) or warn $!;

# print the reversed decoy protein entries
my $current_header = '';
my $current_sequence = '';
open(I, $indb) or die $!;
while(<I>) {
	my $line = $_;
	chomp($line);
	if ( $line =~ /^>/ ) {
		print_reverse( $current_header, $current_sequence, $outfh );
		$current_header = $line;
		$current_sequence = '';
	}
	elsif ( $line =~ /^[A-Z]+$/ ) {
		$current_sequence .= $line;
	}
}
print_reverse( $current_header, $current_sequence, $outfh );

# close the filehandles
close(I) or warn $!;
$outfh->close();


#-------------------------------------------------------------
# SUBS
#-------------------------------------------------------------
sub print_reverse {
	my ( $current_header, $current_sequence, $outfh ) = @_;
	
	my $chpl = 60;
	
	# return if header or sequence is not yet defined
	return if $current_header eq '' or $current_sequence eq '';
	
	my $new_header = change_header( $current_header );
	my $rev_sequence = reverse( $current_sequence );
	my $length = length( $rev_sequence );
	
	# calculate the number of rows needed for the sequence
	my $mod = $length % $chpl;
	my $rows = ( $length - $mod ) / $chpl;
	$rows++ if $mod > 0;
	
	# print the header
	print $outfh $new_header . "\n";
	
	# print the sequence 
	for ( my $i = 0; $i < $rows; $i++ ) {
		my $start = $i * $chpl;
		
		# substr takes the string until the end if $len is too big
		print $outfh substr( $rev_sequence, $start, $chpl ) . "\n";
	}
}

sub change_header {
	my ( $h ) = @_;
	my $new = '';
	
	$h =~ s/^>//;
	$new = '>' . $prefix . $h;
	
	return $new;
}


__END__


Lukas Reiter
