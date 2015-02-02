package BinnedEntity;
use strict;

##################################################################
#
# BinnedEntity
#
# NOTES:
# Object representing binned higher level constructs like
# peptides or proteins
#
# SYNOPSIS:
#
# TODO:
# -
# -
#
# CHANGELOG:
# 2008.02.20
# -
# -
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $at = AnalysisType->new( $verbose, $status );
# Function  :
# Returns   :
# Args      :
sub new {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my ( $verbose, $status ) = @_;

	$self->{ids} = {};
	
	if ( defined($verbose) ) {
		$self->{v} = $verbose;
	}
	else {
		$self->{v} = 0;
	}
	if ( defined($status) ) {
		$self->{s} = $status;
	}
	else {
		$self->{s} = 0;
	}

	return $self;
}

# Object Method
# Title	    :  add_feature()
# Usage     :
# Function	:  ads a feature (e.g. bin) value pair
#              for an id (e.g. protein id)
# Returns   :
# Args      :
sub add_feature {
	my $self = shift;
	my ( $ra_ra_id_feature, $feature_name ) = @_;

	# {ids}{entity_id}{feature_name(e.g. bin)} = value
	foreach (@$ra_ra_id_feature) {
		$self->{ids}{ $_->[0] }{$feature_name} = $_->[1];
	}
}

# Object Method
# Title	    :  add_feature()
# Usage     :
# Function	:  ads a feature (e.g. bin) value pair
#              for an id (e.g. protein id)
# Returns   :
# Args      :
sub get_feature_by_id {
	my $self = shift;
	my ( $id, $feature_name ) = @_;

	if ( exists($self->{ids}{$id}{$feature_name}) ) {
		return $self->{ids}{$id}{$feature_name};
	}
	else {
		print "  feature $feature_name or id $id not found!\n";
		return 0;
	}
}

# Object Method
# Title	    :  add_bins_ob()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub add_bins_ob {
	my $self = shift;
	my ($bins_ob) = @_;

	$self->{bins} = $bins_ob;
}

# Object Method
# Title	    :  get_nr_bins()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nr_bins {
	my $self = shift;
	
	if ( exists($self->{bins}) ) {
		return $self->{bins}->get_nr_bins();
	}
	else {
		print "  no bins found!\n";
		return 0;
	}
}

# Object Method
# Title	    :  get_bin_by_index()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_bin_by_index {
	my $self = shift;
	my ($index) = @_;

	if ( exists($self->{bins}) ) {
		return $self->{bins}->get_bin_by_index($index);
	}
	else {
		print "  no bins found!\n";
		return 0;
	}
}

# Object Method
# Title	    :  get_nr_ids()
# Usage     :
# Function	:  returns the number of ids that were found
#              in bin $i
# Returns   :
# Args      :
sub get_nr_ids {
	my $self = shift;
	my ( $ra_ids, $i ) = @_;

	my $nr_ids = 0;
	foreach (@$ra_ids) {
		if ( exists( $self->{ids}{$_}{bin} ) ) {
			if ( $self->{ids}{$_}{bin} == $i ) {
				$nr_ids++;
			}
		}
	}

	return $nr_ids;
}

# Object Method
# Title	    :  get_ids()
# Usage     :
# Function	:  
# Returns   :
# Args      :
sub get_ids {
	my $self = shift;
	
	my $rh_ids;
	
	foreach my $id ( keys% {$self->{ids}} ) {
		$rh_ids->{$id}++;
	}

	return $rh_ids;
}

# Object Method
# Title	    :  id_exists()
# Usage     :
# Function	:  
# Returns   :
# Args      :
sub id_exists {
	my $self = shift;
	my ( $id ) = @_;
	
	if ( exists($self->{ids}{$id})) {
		return 1;
	}
	else {
		return 0;
	}
}

sub n_print {
	print "  " . $_[0];
}

sub v_print {
	my $self = shift();
	print "  " . $_[0] if $self->{v};
}

1;

__END__

Lukas Reiter

