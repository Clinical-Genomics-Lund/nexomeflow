#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use List::Util qw( min max );
my %vcf_meta;
my @vcf_data;
my @head;
open(VEP, $ARGV[0]);
my $c = 0;
while( <VEP>) {
    ## Print and store Meta-info
    if( /^##/ ) {
        print;
        $c++;
        if ( $c == 10 ) {
        print "##INFO\=<ID=GNOMADAF\,Number=1\,Type=Float,Description=\"Average AF GnomAD\">\n";
        print "##INFO=<ID=GNOMADAF_MAX,Number=1,Type=Float,Description=\"Highest reported AF in gnomAD\">\n";
        print "##INFO=<ID=dbNSFP_GERP___RS,Number=1,Type=Float,Description=\"GERP score\">\n";
        print "##INFO=<ID=dbNSFP_phyloP100way_vertebrate,Number=1,Type=Float,Description=\"phyloP100 score\">\n";
        print "##INFO=<ID=dbNSFP_phastCons100way_vertebrate,Number=1,Type=Float,Description=\"phastcons score\">\n";
        }
        my( $type, $meta ) = parse_metainfo( $_ );
	    $vcf_meta{$type}->{$meta->{ID}} = $meta if defined $type;


    }
    # Print and store header
    elsif( /^#/ ) {
	    print;
        $_ =~ s/^#//;
	    @head = split /\t/;

    }
    # Print and store variant information, add gnomadg and conservation scores
    # to info-field.
    else {
        my $doobi = parse_variant( $_, \@head, \%vcf_meta );
        my @add_info_field;

        my @VARIANTS = split /\t/;
        print join "\t", @VARIANTS[0..6];

        my @info_field = split/;/,$VARIANTS[7];
        print "\t";

        my $gAF = $doobi->{INFO}->{CSQ}->[0]->{gnomADg_AF};
        if ($gAF) {
             push @add_info_field,"GNOMADAF=$gAF";
        }
        my @sub_pop_af;
        foreach my $subpop (keys $doobi->{INFO}->{CSQ}->[0]) {
            if ($subpop =~ /gnomADg_AF_/) {
                push  @sub_pop_af,$doobi->{INFO}->{CSQ}->[0]{$subpop};
            }
        }
        my $max = findmax(@sub_pop_af);
        if ($max > 0) {
            push @add_info_field,"GNOMADAF_MAX=$max";
        }

        my $GERP = $doobi->{INFO}->{CSQ}->[0]->{GERP};
        if ($GERP) {
            push @add_info_field,"dbNSFP_GERP___RS=$GERP";
        }
        my $pC = $doobi->{INFO}->{CSQ}->[0]->{phastCons};
        if ($pC) {
            push @add_info_field,"dbNSFP_phastCons100way_vertebrate=$pC";

        }
        my $pP = $doobi->{INFO}->{CSQ}->[0]->{phyloP100way};
         if ($pP) {
            push @add_info_field,"dbNSFP_phyloP100way_vertebrate=$pP";
        }
        #Add new info field information
        push @info_field, @add_info_field;
        #print new and old information
        print join ";", @info_field;
        print "\t";
        #print everything after info field
        print join "\t", @VARIANTS[8..$#VARIANTS];
    }
}

sub findmax {
    my @in = @_;
    my $high = 0;
    foreach my $val (@in) {
        if ($val && $val ne '.') {
           if($val > $high) {
               $high = $val;
           }
       }
       else {

       }
    }
    return $high;
}

# Parse VCF meta info line (only FORMAT and INFO)
sub parse_metainfo {
    my $comment = shift;

    $comment =~ s/^##//;
    my( $type, $data ) = ( $comment =~ /^(.*?)=(.*)$/ );


    if( $type eq "FORMAT" or $type eq "INFO" or $type eq "SAMPLE" or $type eq "FILTER" ) {
	$data = remove_surrounding( $data, '<', '>' );
	my $pairs = keyval( $data, '=', ',' );
	return $type, $pairs;
    }

    return undef, undef;
}


# Parse VCF variant line
sub parse_variant {
    my( $var_str, $head, $meta ) = @_;
    my @var_data = split /\t/, $var_str;
    my %var;

    $var{ vcf_str } = $var_str;

    # First seven fields
    for ( 0..6 ) {
	$var{ $head->[$_] } = $var_data[$_];
    }

    # Eigth field, INFO
    $var{ INFO } = parse_info( $var_data[7] );

    # Parse VEP annotation field, if any
    if( $var{ INFO }->{ CSQ } ) {
	$var{ INFO }->{ CSQ } = parse_VEP_CSQ( $var{INFO}->{CSQ}, $meta->{INFO}->{CSQ} );
    }

    # Genotypes for each sample
    for ( 9 .. (@var_data-1) ) {
	$var{ GT } -> { $head->[$_] } = parse_genotype( $var_data[8], $var_data[$_] );
    }

    return \%var;
}


# Parse genotype field of VCF
sub parse_genotype {
    my( $format, $data ) = @_;

    my @format = split ':', $format;
    my @data   = split ':', $data;

    my %gt;
    @gt{@format} = @data;

    return \%gt;
}


# Parse info column of VCF file
sub parse_info {
    my $str = shift;
    my $info = keyval( $str, "=", ";" );

    return $info;
}


sub parse_VEP_CSQ {
    my( $CSQ_var, $CSQ_meta ) = @_;

    $CSQ_meta->{Description} =~ /Consequence annotations from Ensembl VEP\. Format: (.*?)$/;

    my @field_names = split /\|/, $1;

    my @transcripts = split /,/, $CSQ_var;

    my @data_transcripts;
    foreach my $transcript_CSQ ( @transcripts ) {
	my @values = split /\|/, $transcript_CSQ;

	my %data;
	for( 0 .. $#field_names ) {
	    if( $field_names[$_] eq "Consequence" ) {
		my @conseq_array = split '&', $values[$_];
		$data{ $field_names[$_] } = \@conseq_array;
	    }
	    else {
		$data{ $field_names[$_] } = ( $values[$_] or "" );
	    }
	}

	push( @data_transcripts, \%data )
    }
    return \@data_transcripts;
}

# Removes character(s) defined in arg2 if first in string, and arg3 if last in string.
sub remove_surrounding {
    my( $str, $before, $after ) = @_;
    $str =~ s/^$before//;
    $str =~ s/$after$//;
    return $str;
}


# Parse string with key value pairs. Return hash.
#  * Keys and values separated by 2nd argument.
#  * Pairs separated by 3rd argument
#  * Handles commas in values if surrounded by double quotes
sub keyval {
    my( $str, $keyval_sep, $pair_sep ) = @_;

    my @pair_str = split /$pair_sep/, $str;
    my %pairs;
    foreach( @pair_str ) {

	# If key-value separator exists, save the value for the key
	if( /$keyval_sep/ ) {
	    my( $key, $val ) = split /$keyval_sep/;
        $val = remove_surrounding( $val, '"', '"' );
	    $pairs{$key} = $val;
	}

	# Otherwise treat the whole string as a flag and set it to one (true).
	else {
	    $pairs{$_} = 1;
	}
    }
    return \%pairs;
}

sub excel_float {
    my $val = shift;

    return 0 if $val eq ".";

    $val =~ s/\./,/;
    return $val;
}
