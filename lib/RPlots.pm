package RPlots;
use strict;

use File::Basename;

##################################################################
#
# RPlots
#
#
# NOTES:
#
# SYNOPSIS:
#
# TODO:
# -
# -
#
# CHANGELOG:
# 2008.04.15
# -
# -
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $r = RPlots->new( $verbose, $status );
# Function  :
# Returns   :
# Args      :
sub new {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my ( $verbose, $status ) = @_;

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

	$self->{r_command} = "R";

	return $self;
}

# Object Method
# Title	    :  run_r_template()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub run_r_template {
	my $self = shift();
	my ( $template_in, $rh_replaces, $template_out_dir ) = @_;

	if ( defined($template_in) ) {
		if ( -e $template_in ) {

			my $r_script_name = basename($template_in);
			if ( defined($template_out_dir) ) {
				if ( -e $template_out_dir ) {
					$r_script_name =
					  $template_out_dir . basename($template_in);
				}
			}

			# read
			open( T, $template_in ) or warn $!;
			my @template = <T>;
			close(T);

			# replace
			my @replace = ();
			foreach my $line (@template) {
				chomp($line);
				foreach my $pattern ( keys %$rh_replaces ) {
					my $substitute = $rh_replaces->{$pattern};
					if ( $line =~ /$pattern/ ) {
						$line =~ s/$pattern/$substitute/;
					}
				}
				push @replace, $line;
			}

			# write
			open( O, ">$r_script_name" ) or warn $!;
			print O join( "\n", @replace ) . "\n";
			close(O);

			# run
			my $command = $self->{r_command} . " CMD BATCH $r_script_name";
			$self->v_print("R command:\n");
			$self->v_print("$command\n");
			system($command );

		}
		else {
			print "  '$template_in' not found!\n";
			return;
		}
	}
	else {
		print "  '$template_in' not defined!\n";
		return;
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





