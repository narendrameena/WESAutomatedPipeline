# author:       narumeena
# description:  devlop automaticaly shell script of variant calling pipeline
# purpose:      read a tab-delimited text file (i.e., the "fields" are separated by
#               the tab \t
# usage:        perl main.pl


#!usr/bin/perl -w

$filename = 'sample.txt';

use strict;
use Cwd;
use Cwd 'abs_path';
use Getopt::Long;
use Config::Abstract::Ini;


###############################################################################
# Command line options
#
my $versionString = "1.0";

my $miniManual = qq"
Usage: perl main.pl [-h] -d <dataTable> -c <config> -o <outputDir>

Options:
[-v|--version]                                 : Print out just the version.
[-h|--help]                                    : Print this miniManual.
[-d|--dataTable]   <tab-delimited file>        : Data table with all samples
[-c|--config]      <configure file>            : Configure file.
[-o|--outputDir]   <output directory>          : Output directory.


";

if ((@ARGV == 0) || ($ARGV[0] eq "-h")||($ARGV[0] eq "--help")) {
    print "$miniManual\n";
    exit();
}

my $dataTable = "";
my $outputDir = "";
my $config = "";
my $command = "";
my $version=0;

GetOptions(
"help|h",
"version|v" => \$version,
"dataTable|d:s" => \$dataTable,
"config|c:s" => \$config,
"outputDir|o:s" => \$outputDir,
) or die "Unknown option\n";

if ($version) {
    print "Version: $versionString\n";
    exit();
}

$outputDir =~ s/\/$//g;



###############################################################################
# Config file
#
my $ini = new Config::Abstract::Ini($config);
my %paths = $ini->get_entry('PATHS');
my %parameters = $ini->get_entry('PARAMETERS');

my $fastqc=$paths{'fastqc'};
my $bowtie2=$paths{'bowtie2'};
my $samtools=$paths{'samtools'};
my $sortsam=$paths{'sortsam'};
my $mergsamfile=$paths{'mergsamfile'};
my $markduplicate=$paths{'markduplicate'};
my $GATK=$paths{'GATK'};
my $hg19Fasta=$paths{'hg19Fasta'};
my $MillsIndels=$paths{'MillsIndels'};
my $hapmap=$paths{'hapmap'};
my $omni=$paths{'omni'};
my $G1000=$paths{'G1000'};
my $ThousandGphaseIndels=$paths{'ThousandGphaseIndels'};
my $dbsnp138=$paths{'dbsnp138'};
my $convert2annovar=$paths{'convert2annovar'};
my $summarizeannovar=$paths{'summarizeannovar'};
my $annovardatabase=$paths{'annovardatabase'};
my $VEP=$paths{'VEP'};
my $bgzip=$paths{'bgzip'};
my $tabix=$paths{'tabix'};
my $cluster=$paths{'cluster'};

my $fastQCFormat=$parameters{'fastQCFormat'};
my $picardMem=$parameters{'picardMem'};
my $gtkMem=$parameters{'gtkMem'};
my $picardAs=$parameters{'picardAs'};
my $picardUsethread=$parameters{'picardUsethread'};
my $markduplicateRemoveduplicate=$parameters{'markduplicateRemoveduplicate'};
my $markduplicateAs=$parameters{'markduplicateAs'};
my $markduplicateValidationstringency=$parameters{'markduplicateValidationstringency'};
my $unifiedGenotyperStandcallconf=$parameters{'unifiedGenotyperStandcallconf'};
my $unifiedGenotyperStandemitconf=$parameters{'unifiedGenotyperStandemitconf'};
my $unifiedGenotyperGlm=$parameters{'unifiedGenotyperGlm'};
my $unifiedGenotyperOutmode=$parameters{'unifiedGenotyperOutmode'};
my $buildver=$parameters{'buildver'};
my $verdbsnp=$parameters{'verdbsnp'};
my $ver1000g=$parameters{'ver1000g'};
my $veresp=$parameters{'veresp'};
my $VQSR=$parameters{'VQSR'};
my $parallel=$parameters{'parallel'};


###############################################################################
# Read datatable
#
open (IN, $dataTable);
#
# Read the title and figure out column information
#
my $title=<IN>;
chomp $title;
$title=~s/\r//g;
my @labels=split(/\s+/,$title);

my %labelsHash=();
for (my $i=0; $i<=$#labels; $i++){
    $labelsHash{$labels[$i]}=$i;
}

# check header columns
my @headerCols = ("Sample", "SeqPlatform", "Library", "SeqType", "ReadGroup", "Directory", "SeqFile1", "SeqFile2");
for(my $i=0; $i<=$#headerCols; $i++){
    if(! exists($labelsHash{$headerCols[$i]})){
        print STDERR "Oops, a column is missing in the data table header: $headerCols[$i].\n";
        exit();
    }
}



###############################################################################
# Read the data
#
my ($Sample, $SeqPlatform, $Library, $SeqType, $ReadGroup, $Directory, $SeqFile1, $SeqFile2);
my $path="";
my %SampleBamMap=();

for my $line (<IN>){
    chomp $line;
    next if $line =~ /^\s*$/; # skip blank lines
    $line =~ s/\r//g;
    my @data = split(/\s+/,$line);
    $Sample=$data[$labelsHash{Sample}];
    $SeqPlatform=$data[$labelsHash{SeqPlatform}];
    $Library=$data[$labelsHash{Library}];
    $SeqType=$data[$labelsHash{SeqType}];
    $ReadGroup=$data[$labelsHash{ReadGroup}];
    $Directory=$data[$labelsHash{Directory}];
    $SeqFile1=$data[$labelsHash{SeqFile1}];
    $SeqFile2=$data[$labelsHash{SeqFile2}];
    $path = "$outputDir\/$Sample";
    print $path,"\n";
    unless (-d $path) {system("mkdir -p $path")};
    
    if (!exists($SampleBamMap{$Sample})){
        $SampleBamMap{$Sample} = "$path\/bowtie2\_$Sample\_bowtie2.bam";
    } else {
        $SampleBamMap{$Sample} .= ",$path\/bowtie2\_$Sample\_bowtie2.bam";
    }
    
    if ($SeqType eq "Paired-end"){
        
        open(QC_Mapping, ">$outputDir\/QC_Mapping\_$Sample\_bowtie2\.sh");
        print QC_Mapping "\#!\/bin\/sh\n\n";
        if ($parallel==1){
            open(Cluster, "<$cluster");
            while (<Cluster>){
                print QC_Mapping $_;
            }
            close Cluster;
        }
        print QC_Mapping "\#Run QC for Sample $Sample \(Read group: $ReadGroup\)\:\n\n";
        print QC_Mapping "$fastqc -o $path -f $fastQCFormat $Directory\/$SeqFile1 $Directory\/$SeqFile2 \>\& $path\/fastqc\_$Sample\_bowtie2\.log\n\n";
        print QC_Mapping "$bowtie2  $hg19Fasta $Directory\/$SeqFile1 $Directory\/$SeqFile2  > $path\/bowtie2\_$Sample\_bowtie2.sam 2\> $path\/bowtie2\_$Sample\_bowtie2.log\n\n";
        print QC_Mapping "$samtools view -bS $path\/bowtie2\_$Sample\_bowtie2.sam >$path\/bowtie2\_$Sample\_bowtie2.bam 2>$path\/bowtie2\_$Sample\_bowtie2\_samtobam.log\n\n";
        close QC_Mapping;
        
    }elsif ($SeqType eq "Single-end") {
        
        open(QC_Mapping, ">$outputDir\/QC_Mapping\_$Sample\_bowtie2\.sh");
        print QC_Mapping "\#!\/bin\/sh\n\n";
        if ($parallel==1){
            open(Cluster, "<$cluster");
            while (<Cluster>){
                print QC_Mapping $_;
            }
            close Cluster;
        }
        print QC_Mapping "\#Run QC for Sample $Sample \(Read group: $ReadGroup\)\:\n\n";
        print QC_Mapping "$fastqc -o $path -f $fastQCFormat $Directory\/$SeqFile1 \>\& $path\/fastqc\_$Sample\_bowtie2\.log\n\n";
        print QC_Mapping "$bowtie2  $hg19Fasta $Directory\/$SeqFile1 > $path\/bowtie2\_$Sample\_bowtie2.sam 2\> $path\/bowtie2\_$Sample\_bowtie2.log\n\n";
        print QC_Mapping "$samtools view -bS $path\/bowtie2\_$Sample\_bowtie2.sam >$path\/bowtie2\_$Sample\_bowtie2.bam 2>$path\/bowtie2\_$Sample\_bowtie2\_samtobam.log\n\n";
        close QC_Mapping;
    }else{
        print "Oops, the SeqType is not set at a correct value!";
    }
}
close IN;

