#!/usr/bin/perl

#   adapted from FASTX-toolkit -fastx_barcode_splitter
#

use strict;
use warnings;
use IO::Handle;
use Data::Dumper;
use Getopt::Long;
use Carp;
#use PerlIO::gzip;

##
## Based on barcode matching.
##
## run with "--help" for usage information
##

# Forward declarations
sub load_barcode_file ($);
sub parse_command_line ;
sub match_sequences ;
sub mismatch_count($$) ;
sub print_results;
sub open_and_detect_input_format;
sub read_record;
sub write_record($);
sub usage();

# Global flags and arguments, 
# Set by command line argumens
my $barcode_file ;
my $sample_file ;
my $barcodes_at_eol = 0 ;
my $barcodes_at_bol = 0 ;
my $exact_match = 0 ;
my $allow_partial_overlap = 0;
my $allowed_mismatches = 1;
my $newfile_suffix = '';
my $newfile_prefix  ;
my $quiet = 0 ;
my $debug = 0 ;
my $fastq_format = 1;

# Global variables 
# Populated by 'create_output_files' for r1 r2
my %filenames1;
my %filenames2;
my %files1;
my %files2;
my %bar_seq;

my %counts = ( 'unmatched' => 0 );
my $barcodes_length;
my @barcodes;
my @sample_info;
my $input_file_io;


# The Four lines per record in FASTQ format.
# (when using FASTA format, only the first two are used)
my $seq_name;
my $seq_bases;
my $seq_name2;
my $seq_qualities;

#get options
my ($indir,$outdir,$sample,$quality,$end);

my $R1="_R1.fq.gz";
my $R2="_R2.fq.gz";

# Start of Program
#
#
parse_command_line ;

load_barcode_file ( $barcode_file ) ;


#match_sequences ;

#print_results unless $quiet;

#
# End of program
#








sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV==0);

	my $result = GetOptions ( "bcfile=s" => \$barcode_file,
				  "infile=s" => \$sample_file,
				  "eol"  => \$barcodes_at_eol,
				  "bol"  => \$barcodes_at_bol,
				  "sample=s"=>\$sample,
				  "indir=s"=>\$indir,
				  "outdir=s"=>\$outdir,
				  "exact" => \$exact_match,
				  "prefix=s" => \$newfile_prefix,
				  "suffix=s" => \$newfile_suffix,
				  "quiet" => \$quiet, 
				  "partial=i" => \$allow_partial_overlap,
				  "debug" => \$debug,
				  "mismatches=i" => \$allowed_mismatches,
				  "help" => \$help,
				  "end=s"=>\$end,

				  ) ;
	
	usage() if ($help);

	die "Error: barcode file not specified (use '--bcfile [FILENAME]')\n" unless defined $barcode_file;
	die "Error: prefix path/filename not specified (use '--prefix [PATH]')\n" unless defined $newfile_prefix;
        die "Error: sample in file not specified (use '--infile [FILENAME]')\n" unless defined $sample_file;
	if ($barcodes_at_bol == $barcodes_at_eol) {
		die "Error: can't specify both --eol & --bol\n" if $barcodes_at_eol;
		die "Error: must specify either --eol or --bol\n" ;
	}

	die "Error: invalid for value partial matches (valid values are 0 or greater)\n" if $allow_partial_overlap<0;

	$allowed_mismatches = 0 if $exact_match;

	die "Error: invalid value for mismatches (valid values are 0 or more)\n" if ($allowed_mismatches<0);

	die "Error: partial overlap value ($allow_partial_overlap) bigger than " . 
		"max. allowed mismatches ($allowed_mismatches)\n" if ($allow_partial_overlap > $allowed_mismatches);


	exit unless $result;
}

$end ||=2;
$indir = trim_slash($indir);
`mkdir -p $newfile_prefix` unless(-d "$newfile_prefix");



#
# Read the barcode file
#
sub load_barcode_file ($) {
	my $filename = shift or croak "Missing barcode file name";

	open BCFILE,"<$filename" or die "Error: failed to open barcode file ($filename)\n";
	while (<BCFILE>) {
		next if m/^#/;
		chomp;
		my ($ident, $barcode) = split ;

		$barcode = uc($barcode);

		# Sanity checks on the barcodes
		die "Error: bad data at barcode file ($filename) line $.\n" unless defined $barcode;
		die "Error: bad barcode value ($barcode) at barcode file ($filename) line $.\n"
			unless $barcode =~ m/^[AGCT]+$/;

		die "Error: bad identifier value ($ident) at barcode file ($filename) line $. (must be alphanumeric)\n" 
			unless $ident =~ m/^\w+$/;

		die "Error: badcode($ident, $barcode) is shorter or equal to maximum number of " .
		    "mismatches ($allowed_mismatches). This makes no sense. Specify fewer  mismatches.\n" 
		    	if length($barcode)<=$allowed_mismatches;

		$barcodes_length = length($barcode) unless defined $barcodes_length;
		die "Error: found barcodes in different lengths. this feature is not supported yet.\n" 
			unless $barcodes_length == length($barcode);

	 	push @barcodes, [$ident, $barcode];
		$bar_seq{$ident} = $barcode;
		if ($allow_partial_overlap>0) {
			foreach my $i (1 .. $allow_partial_overlap) {
				substr $barcode, ($barcodes_at_bol)?0:-1, 1, '';
	 			push @barcodes, [$ident, $barcode];
				$bar_seq{$ident} = $barcode;
			}
		}
	}
	close BCFILE;

	if ($debug) {
		print STDERR "barcode\tsequence\n";
		foreach my $barcoderef (@barcodes) {
			my ($ident, $seq) = @{$barcoderef};
			print STDERR $ident,"\t", $seq ,"\n";
		}
	}
}


# Create one output file for each barcode.
# (Also create a file for the dummy 'unmatched' barcode)
sub create_output_files {
	my %barcodes = map { $_->[0] => 1 } @sample_info; #generate a uniq list of barcode identifiers(revised to gennerate files by sample list provided);
	$barcodes{'unmatched'} = 1 ;

	foreach my $ident (keys %barcodes) {
	#	my $sample_name =$samplenames
		my $new_filename1 = $newfile_prefix . $ident . $newfile_suffix . $R1;
		my $new_filename2 = $newfile_prefix . $ident . $newfile_suffix . $R2;
		$filenames1{$ident} = $new_filename1;
		$filenames2{$ident} = $new_filename2;
		open my $file1, "|gzip > $new_filename1" or die "Error: failed to create output file ($new_filename1)\n"; 
         	open my $file2, "|gzip > $new_filename2" or die "Error: failed to create output file ($new_filename2)\n";
		$files1{$ident} = $file1 ;
		$files2{$ident} = $file2 ;
	}
}
#------subroutine to get sample info ------
open FILE,"<$sample_file" or die "Error: failed to open sample file ($sample_file)\n";
while (<FILE>) {
	next if m/^#/;
	chomp;
	my @f = split /\s+/;
	my ($sample_name, $ident)=($f[2],$f[1]);
	push @sample_info, [$sample_name, $ident];	
	#push @barcodes, [$ident, $barcode];
}


#-----------------main project to output 
if ($end == 2) {

        #-get the input files

        chomp (my $file_1 = `ls $indir/${sample}/*_*1.fq.gz`);
        chomp (my $file_2 = `ls $indir/${sample}/*_*2.fq.gz`);

        #-open the input files

        open IN_1,"gzip -dc $file_1 |" or die $!;
        #open IN_2,"<:gzip","$file_2" or die $!;
        open IN_2,"gzip -dc $file_2 |" or die $!;
	create_output_files;
       
        #---------------------------read in the reads information-----------------------#

        while (1) {
                #-get the reads and corresponding information in each 4 lines

                my $line1_1 = <IN_1>;
                my $line1_2 = <IN_1>;
                my $line1_3 = <IN_1>;
                my $line1_4 = <IN_1>;

                my $line2_1 = <IN_2>;
                my $line2_2 = <IN_2>;
                my $line2_3 = <IN_2>;
                my $line2_4 = <IN_2>;

                #check the end of the file

                last unless (defined($line1_1) and defined($line2_1));
                chomp ($line1_1,$line1_2,$line1_3,$line1_4,$line2_1,$line2_2,$line2_3,$line2_4);
                my $best_barcode_mismatches_count = $barcodes_length;
                my $best_barcode_ident = undef;
                #Try all barcodes, find the one with the lowest mismatch count
		foreach my $sample (@sample_info){
		        my ($sample_id,$ident) =@{$sample};
			my $barcode = $bar_seq{$ident};
                       # my ($ident, $barcode) = @{$barcoderef};

                        # Get DNA fragment (in the length of the barcodes)
                        # The barcode will be tested only against this fragment
                        # (no point in testing the barcode against the whole sequence)
                        my $sequence_fragment;
                        if ($barcodes_at_bol) {
                                $sequence_fragment = substr $line2_2, 0, $barcodes_length;
                        } else {
                                $sequence_fragment = substr $line2_2, - $barcodes_length;
                        }
                        my $mm = mismatch_count($sequence_fragment, $barcode) ;

                        # if this is a partial match, add the non-overlap as a mismatch
                        # (partial barcodes are shorter than the length of the original barcodes)
                        $mm += ($barcodes_length - length($barcode));

                        if ( $mm < $best_barcode_mismatches_count ) {
                                $best_barcode_mismatches_count = $mm ;
                                $best_barcode_ident = $sample_id ;
                        }
                }

                $best_barcode_ident = 'unmatched'
                        if ( (!defined $best_barcode_ident) || $best_barcode_mismatches_count>$allowed_mismatches) ;

                print STDERR "sequence $seq_bases matched barcode: $best_barcode_ident\n" if $debug;

                $counts{$best_barcode_ident}++;
                my $file1 = $files1{$best_barcode_ident};
                my $file2 = $files2{$best_barcode_ident};
	     	if ( $best_barcode_ident ne 'unmatched') {
		#	print "$best_barcode_ident\n";
			print $file1 "$line1_1\n$line1_2\n$line1_3\n$line1_4\n";
               		print $file2 "$line2_1\n$line2_2\n$line2_3\n$line2_4\n";
		}

        }
}


print_results;


#Quickly calculate hamming distance between two strings
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }



sub print_results
{
	my $log = "_barcode_number.log";
	my $barcode_filename = $newfile_prefix . $sample . $log;
	open my $file1, " > $barcode_filename" or die "Error: failed to create output file ($barcode_filename)\n";
	print $file1 "Barcode\tCount\tLocation\n";
	my $total = 0 ;
	foreach my $ident (sort keys %counts) {
		print $file1 "$ident\t$counts{$ident}\t$filenames1{$ident}\n";
		$total += $counts{$ident};
	}
	print  $file1 "total\t$total\n";
}


sub read_record
{
	$seq_name = $input_file_io->getline();

	return undef unless defined $seq_name; # End of file?

	$seq_bases = $input_file_io->getline();
	die "Error: bad input file, expecting line with sequences\n" unless defined $seq_bases;

	# If using FASTQ format, read two more lines
	if ($fastq_format) {
		$seq_name2  = $input_file_io->getline();
		die "Error: bad input file, expecting line with sequence name2\n" unless defined $seq_name2;

		$seq_qualities = $input_file_io->getline();
		die "Error: bad input file, expecting line with quality scores\n" unless defined $seq_qualities;
	}
	return 1;
}

sub write_record($)
{
	my $file = shift;

	croak "Bad file handle" unless defined $file;

	print $file $seq_name;
	print $file $seq_bases,"\n";

	#if using FASTQ format, write two more lines
	if ($fastq_format) {
		print $file $seq_name2;
		print $file $seq_qualities;
	}
}

#open file
sub open_and_detect_input_format
{
	$input_file_io  = new IO::Handle;
	die "Failed to open STDIN " unless $input_file_io->fdopen(fileno(STDIN),"r");

	# Get the first characeter, and push it back
	my $first_char = $input_file_io->getc();
	$input_file_io->ungetc(ord $first_char);

	if ($first_char eq '>') {
		# FASTA format
		$fastq_format = 0 ;
		print STDERR "Detected FASTA format\n" if $debug;
	} elsif ($first_char eq '@') {
		# FASTQ format
		$fastq_format = 1;
		print STDERR "Detected FASTQ format\n" if $debug;
	} else {
		die "Error: unknown file format. First character = '$first_char' (expecting > or \@)\n";
	}
}


#-dir trimming
sub trim_slash {
	my($dir) = @_;
	($dir =~ /\/$/) ? ($dir =~ s/\/$//) : ($dir = $dir);
        return $dir;
                }
                        



sub usage()
{
print<<EOF;

usage: $0 --bcfile FILE --prefix PREFIX [--suffix SUFFIX] [--bol|--eol] 
         [--mismatches N] [--exact] [--partial N] [--help] [--quiet] [--debug]

Arguments:

--bcfile FILE	- Barcodes file name. (see explanation below.)

--infile FILE   - SAMPLE list file .

--prefix PREFIX	- File prefix. will be added to the output files. Can be used
		  to specify output directories.
--suffix SUFFIX	- File suffix (optional). Can be used to specify file
		  extensions.
--bol		- Try to match barcodes at the BEGINNING of sequences.
		  (What biologists would call the 5' end, and programmers
		  would call index 0.)
--eol		- Try to match barcodes at the END of sequences.
		  (What biologists would call the 3' end, and programmers
		  would call the end of the string.)

		  NOTE: one of --bol, --eol must be specified, but not both.
--mismatches N	- Max. number of mismatches allowed. default is 1.
--exact		- Same as '--mismatches 0'. If both --exact and --mismatches 
		  are specified, '--exact' takes precedence.
--partial N	- Allow partial overlap of barcodes. (see explanation below.)
		  (Default is not partial matching)
--quiet		- Don't print counts and summary at the end of the run.
		  (Default is to print.)
--debug		- Print lots of useless debug information to STDERR.
--help		- This helpful help screen.


Barcode file format
-------------------

##   Sample Barcode Rename  ERCC
#    1       1       F_BL_1  10000
     2       2       F_BL_2  10000
     3       3       F_BL_3  10000
     4       4       F_BL_4  10000

EOF

exit 1;
}
