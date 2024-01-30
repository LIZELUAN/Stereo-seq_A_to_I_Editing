#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# A-to-I RNA editing calling based on BAM created using HISAT2
# Jinrong Huang created it for A-to-I RNA editing of bulk RNA-seq data
# Zeluan Li revised it to make it work in stereo-seq data
# Nov 2022

my ($dataset,$sam,$outdir,$samtools);
my ($phred,$qual_cutoff,$depth);
my $suffix;

GetOptions(
	"dataset=s" => \$dataset,
	"sam=s" => \$sam,
	"suffix=s" => \$suffix,
	"outdir=s" => \$outdir,
	"samtools=s" => \$samtools,
	"phred=i" => \$phred,
	"qual_cutoff=i" => \$qual_cutoff,
	"depth=i" => \$depth,
);

$samtools ||="samtools";
$phred ||="33";
$qual_cutoff ||="30";
$depth ||="4";
$suffix ||="sam";

`mkdir -p $outdir` unless (-e $outdir);

# A-to-G
my %hash = (
	"A" => "G",
	"T" => "C",
);

my %dataset;	# known RNA editing sites
if($dataset=~/\.gz$/){open IN,"gunzip -cd <$dataset|" or die $!;}else{open IN,"<$dataset" or die $!;}
while (<IN>){
	# Chromosome_ID   Coordinate      Ref_base        CodonChange     AminoAcidChange Gene    Annotation      Detailed information    Repeat
	# 1       25274   A       -       -       ENSSSCG00000027257:PSMB1        intronic        ENSSSCG00000027257:PSMB1:intronic       SINE/tRNA
	# 1       4846284 T       -       -       ENSSSCG00000004029:QKI  intronic        ENSSSCG00000004029:QKI:intronic SINE/tRNA
#	next if /^Chromosome/;
	chomp;
	next if $.==1;
	my ($chr,$pos,$ref)=(split /\t/)[0,1,2];
	$dataset{"$chr\t$pos"}=$ref;
}
close IN;

my $name=(split /\//,$sam)[-1];
$name=~s/\.$suffix$//;

my %sites;
open IN,"<$sam" or die $!;
while (<IN>){
	chomp;
	my ($CX,$CY,$CHR,$POS,$CIGAR,$SEQ,$QUAL)=(split /\t/)[0,1,4,5,7,11,12];
	my @SEQ =split //,$SEQ;
	my @QUAL =split //,$QUAL;
	my @CIGAR=$CIGAR=~/(\d+[A-Z])/g;

	my $refPos=$POS;
	my $seqPos =0;
	for(my $i=0;$i<@CIGAR;$i++){
		if($CIGAR[$i]=~/M/){
			$CIGAR[$i]=~s/M//;
			for (my $j=0;$j<$CIGAR[$i];$j++){
				my $currentPos=$refPos+$j;
				my $ind="$CX\t$CY\t$CHR\t$currentPos";
				my $index="$CHR\t$currentPos";
				if (exists $dataset{$index}){
					$sites{$ind}{'SEQ'} .= $SEQ[$seqPos];
					$sites{$ind}{'QUAL'} .= $QUAL[$seqPos];
				}
				$seqPos++;
			}
			$refPos+=$CIGAR[$i];
		}
		elsif ($CIGAR[$i]=~/D/){
			$CIGAR[$i]=~s/D//;
			$refPos+=$CIGAR[$i];
		}
		elsif ($CIGAR[$i]=~/I/){
			$CIGAR[$i]=~s/I//;
			$seqPos+=$CIGAR[$i];
		}
		elsif ($CIGAR[$i]=~/N/){
			$CIGAR[$i]=~s/N//;
			$refPos+=$CIGAR[$i];
		}
		elsif ($CIGAR[$i]=~/S/){
			$CIGAR[$i]=~s/S//;
			$seqPos+=$CIGAR[$i];
		}
		else {die "Incorrect format of CIGAR. Make sure the input is sam format!\n"}
	}
}
close IN;

my $sam2base="$outdir/$name.sam2base.gz";
if($sam2base=~/\.gz$/){open OT,"|gzip >$sam2base" or die $!;}else{open OT,">$sam2base" or die $!;}
#while (my($k,$v)=each %sites){
foreach my $k(sort keys %sites){
	my $seq=$sites{$k}{'SEQ'};
	my $qual=$sites{$k}{'QUAL'};
	my $reads=length($seq);
	my ($chr,$pos)=(split /\t/,$k)[2,3];
	my $in="$chr\t$pos";
	print OT "$k\t$dataset{$in}\t$seq\t$qual\t$reads\n";
}	
close OT;


my ($cov_sum,$alt_sum)=("0","0");
my $input=$sam2base;
my $output="$outdir/$name.RE.gz";
if($output=~/\.gz$/){open OT,"|gzip >$output" or die $!;}else{open OT,">$output" or die $!;}
if($input=~/\.gz$/){open IN,"gunzip -cd <$input|" or die $!;}else{open IN,"<$input" or die $!;}
while (<IN>){
	# 100(Cx)	200(Cy)	1       21484011        T       TTTTT   FJJJJ   5
	# 100(Cx)	200(Cy)	4       21558506        T       TTTTTTTTTTTTTT  F)JJ<FFJJJJJJJ  14
	chomp;
	my ($cx,$cy,$chr,$pos,$refbase,$seq,$qual,$reads)=split /\t/;
	my @base=split //,$seq;
	my @qual=split //,$qual;
	my ($base_new,$qual_new);
	my $cov=0;
	my ($ref,$alt)=("0","0");
	my $other=0;
	for (my $i=0;$i<@base;$i++){
		my $score=ord($qual[$i])-$phred;
		next if $base[$i] eq "N";
		next if $score < $qual_cutoff;
		if ($base[$i] =~/^$refbase$/i){
			$ref++;
			$cov++;
			$base_new.=$base[$i];
			$qual_new.=$qual[$i];
		}
		elsif ($base[$i] =~/^$hash{$refbase}$/i){
			$alt++;
			$cov++;
			$base_new.=$base[$i];
			$qual_new.=$qual[$i];
		}
		else {
			$other++;
		}
	}
	next if $cov < $depth;
	print OT "$cx\t$cy\t$chr\t$pos\t$refbase\t$cov\t$alt\n";
}
close IN;
close OT;

exit;
