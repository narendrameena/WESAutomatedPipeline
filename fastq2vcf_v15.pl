#!/usr/bin/perl -w

###############################################################################
# License information
# 
# Commercial users, please contact < ray.x.gao at gmail.com >
#
# Academic users, free of charge. 
#
# Copyright: Gao Lab 
#
###############################################################################

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
Usage: perl fastq2vcf.pl [-h] -d <dataTable> -c <config> -o <outputDir> 

Options:
   [-v|--version]                                 : Print out just the version.
   [-h|--help]                                    : Print this miniManual.
   [-d|--dataTable]   <tab-delimited file>        : Data table with all samples 
   [-c|--config]      <configure file>            : Configure file.
   [-o|--outputDir]   <output directory>          : Output directory.
	
Version: $versionString	
Release date: 05/15/2014
Copyright: Gao Lab
For questions and comments, please contact Gao Lab.
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
my $bwa=$paths{'bwa'};
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
	    $SampleBamMap{$Sample} = "$path\/bwa\_$Sample\_$ReadGroup.bam";
	} else {
	    $SampleBamMap{$Sample} .= ",$path\/bwa\_$Sample\_$ReadGroup.bam";
	}
	
	if ($SeqType eq "Paired-end"){
	    
	    open(QC_Mapping, ">$outputDir\/QC_Mapping\_$Sample\_$ReadGroup\.sh");
	    print QC_Mapping "\#!\/bin\/sh\n\n";
	    if ($parallel==1){
		open(Cluster, "<$cluster");
		while (<Cluster>){
		    print QC_Mapping $_;
		}
		close Cluster; 
	    }
	    print QC_Mapping "\#Run QC for Sample $Sample \(Read group: $ReadGroup\)\:\n\n";
	    print QC_Mapping "$fastqc -o $path -f $fastQCFormat $Directory\/$SeqFile1 $Directory\/$SeqFile2 \>\& $path\/fastqc\_$Sample\_$ReadGroup\.log\n\n";
	    print QC_Mapping "$bwa mem -R '\@RG\\tID:$ReadGroup\\tPL:$SeqPlatform\\tLB:$Library\\tSM:$Sample' $hg19Fasta $Directory\/$SeqFile1 $Directory\/$SeqFile2  > $path\/bwa\_$Sample\_$ReadGroup.sam 2\> $path\/bwa\_$Sample\_$ReadGroup.log\n\n";
	    print QC_Mapping "$samtools view -bS $path\/bwa\_$Sample\_$ReadGroup.sam >$path\/bwa\_$Sample\_$ReadGroup.bam 2>$path\/bwa\_$Sample\_$ReadGroup\_samtobam.log\n\n";
	    close QC_Mapping;
	    
	}elsif ($SeqType eq "Single-end") {
	    
	    open(QC_Mapping, ">$outputDir\/QC_Mapping\_$Sample\_$ReadGroup\.sh");
	    print QC_Mapping "\#!\/bin\/sh\n\n";
	    if ($parallel==1){
                open(Cluster, "<$cluster");
                while (<Cluster>){
		    print QC_Mapping $_;
		}
		close Cluster;
	    }
	    print QC_Mapping "\#Run QC for Sample $Sample \(Read group: $ReadGroup\)\:\n\n";
	    print QC_Mapping "$fastqc -o $path -f $fastQCFormat $Directory\/$SeqFile1 \>\& $path\/fastqc\_$Sample\_$ReadGroup\.log\n\n";
	    print QC_Mapping "$bwa mem -R '\@RG\\tID:$ReadGroup\\tPL:$SeqPlatform\\tLB:$Library\\tSM:$Sample' $hg19Fasta $Directory\/$SeqFile1 > $path\/bwa\_$Sample\_$ReadGroup.sam 2\> $path\/bwa\_$Sample\_$ReadGroup.log\n\n";
	    print QC_Mapping "$samtools view -bS $path\/bwa\_$Sample\_$ReadGroup.sam >$path\/bwa\_$Sample\_$ReadGroup.bam 2>$path\/bwa\_$Sample\_$ReadGroup\_samtobam.log\n\n";
	    close QC_Mapping;
	}else{
	    print "Oops, the SeqType is not set at a correct value!";
	}
}
close IN;

for my $key (keys %SampleBamMap){
    if ($SampleBamMap{$key}=~/\,/){
	my @bamfiles=split("\,",$SampleBamMap{$key});
	my $command="java $picardMem -jar $mergsamfile";
	foreach my $bamfile(@bamfiles){
	    $command .= " I=$bamfile";
        }
	
	open(GATK, ">$outputDir\/Pre-calling\_$key\.sh");
	print GATK "\#!\/bin\/sh\n\n\n";
	if ($parallel==1){
	    open(Cluster, "<$cluster");
	    while (<Cluster>){
                print GATK $_;
	    }
	    close Cluster;
	}   
	print GATK "$command O=$outputDir\/$key\/$key\_merged.bam AS=$picardAs USE_THREADING=$picardUsethread TMP_DIR=$outputDir\/$key\/tmp >& $outputDir\/$key\/$key\_picardmergesamfiles.log\n\n";
	print GATK "java $picardMem -jar $sortsam I=$outputDir\/$key\/$key\_merged.bam O=$outputDir\/$key\/$key\_merged.sorted.bam SO=coordinate TMP_DIR=$outputDir\/$key\/tmp \>\& $outputDir\/$key\/$key\_merged.sorted.log\n\n";
	print GATK "$samtools index $outputDir\/$key\/$key\_merged.sorted.bam \>\&$outputDir\/$key\/$key\_merged.sorted.bam.index.log\n\n";
	print GATK "java $gtkMem -jar $markduplicate I=$outputDir\/$key\/$key\_merged.sorted.bam O=$outputDir\/$key\/$key\_merged.sorted.markdup.bam REMOVE_DUPLICATES=$markduplicateRemoveduplicate AS=$markduplicateAs METRICS_FILE=$outputDir\/$key\/$key\_merged.sorted.markdup.metrics VALIDATION_STRINGENCY=$markduplicateValidationstringency \>\& $outputDir\/$key\/$key\_merged.sorted.markeddup.log\n\n";
	print GATK "$samtools index $outputDir\/$key\/$key\_merged.sorted.markdup.bam \>\&$outputDir\/$key\/$key\_merged.sorted.markdup.index.log\n\n";
	print GATK "java $gtkMem -jar $GATK -R $hg19Fasta -T RealignerTargetCreator -I $outputDir\/$key\/$key\_merged.sorted.markdup.bam -nt 1 -known $MillsIndels -known $ThousandGphaseIndels -o $outputDir\/$key\/$key\_merged.sorted.markdup.bam.RealignerTargetCreator.list -log $outputDir\/$key\/$key\_merged.sorted.markdup.bam.RealignerTargetCreator.log\n\n";
	print GATK "java $gtkMem -jar $GATK -R $hg19Fasta -T IndelRealigner -I $outputDir\/$key\/$key\_merged.sorted.markdup.bam -targetIntervals $outputDir\/$key\/$key\_merged.sorted.markdup.bam.RealignerTargetCreator.list -known $MillsIndels -known $ThousandGphaseIndels -log $outputDir\/$key\/$key\_merged.sorted.markdup.bam.indelrealigner.log -o $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.bam\n\n";
	print GATK "java $gtkMem -jar $GATK -R $hg19Fasta -T BaseRecalibrator -log $outputDir\/$key\/$key\_merged.sorted.markdup.bam.indelrealigner.baserecalibrator.log -knownSites $dbsnp138 -knownSites $MillsIndels -knownSites $ThousandGphaseIndels -I $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.bam -o $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.table\n\n";
	print GATK "java $gtkMem -jar $GATK -R $hg19Fasta -T PrintReads -log $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.table.printreads.log -I $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.bam -BQSR $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.table -o $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.bam\n\n";
	# print GATK "java $gtkMem -jar $GATK -R $hg19Fasta -T ReduceReads -log $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.table.printreads.datacompression.log -I $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.bam -o $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.reduced.bam\n\n";
	close GATK;
	
    } else {
        open(GATK, ">$outputDir\/Pre-calling\_$key\.sh");
        print GATK "\#!\/bin\/sh\n\n\n";
        
	if ($parallel==1){
	    open(Cluster, "<$cluster");
	    while (<Cluster>){
                print GATK $_;
	    }
	    close Cluster;
        }
	print GATK "cp $SampleBamMap{$key} $outputDir\/$key\/$key\_merged.bam\n\n";
        print GATK "java $picardMem -jar $sortsam I=$outputDir\/$key\/$key\_merged.bam O=$outputDir\/$key\/$key\_merged.sorted.bam SO=coordinate TMP_DIR=$outputDir\/$key\/tmp \>\& $outputDir\/$key\/$key\_merged.sorted.log\n\n";
        print GATK "$samtools index $outputDir\/$key\/$key\_merged.sorted.bam \>\&$outputDir\/$key\/$key\_merged.sorted.bam.index.log\n\n";
        print GATK "java $gtkMem -jar $markduplicate I=$outputDir\/$key\/$key\_merged.sorted.bam O=$outputDir\/$key\/$key\_merged.sorted.markdup.bam REMOVE_DUPLICATES=$markduplicateRemoveduplicate AS=$markduplicateAs METRICS_FILE=$outputDir\/$key\/$key\_merged.sorted.markdup.metrics VALIDATION_STRINGENCY=$markduplicateValidationstringency \>\& $outputDir\/$key\/$key\_merged.sorted.markeddup.log\n\n";
        print GATK "$samtools index $outputDir\/$key\/$key\_merged.sorted.markdup.bam \>\&$outputDir\/$key\/$key\_merged.sorted.markdup.index.log\n\n";
        print GATK "java $gtkMem -jar $GATK -R $hg19Fasta -T RealignerTargetCreator -I $outputDir\/$key\/$key\_merged.sorted.markdup.bam -nt 1 -known $MillsIndels -known $ThousandGphaseIndels -o $outputDir\/$key\/$key\_merged.sorted.markdup.bam.RealignerTargetCreator.list -log $outputDir\/$key\/$key\_merged.sorted.markdup.bam.RealignerTargetCreator.log\n\n";
        print GATK "java $gtkMem -jar $GATK -R $hg19Fasta -T IndelRealigner -I $outputDir\/$key\/$key\_merged.sorted.markdup.bam -targetIntervals $outputDir\/$key\/$key\_merged.sorted.markdup.bam.RealignerTargetCreator.list -known $MillsIndels -known $ThousandGphaseIndels -log $outputDir\/$key\/$key\_merged.sorted.markdup.bam.indelrealigner.log -o $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.bam\n\n";
        print GATK "java $gtkMem -jar $GATK -R $hg19Fasta -T BaseRecalibrator -log $outputDir\/$key\/$key\_merged.sorted.markdup.bam.indelrealigner.baserecalibrator.log -knownSites $dbsnp138 -knownSites $MillsIndels -knownSites $ThousandGphaseIndels -I $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.bam -o $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.table\n\n";
        print GATK "java $gtkMem -jar $GATK -R $hg19Fasta -T PrintReads -log $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.table.printreads.log -I $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.bam -BQSR $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.table -o $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.bam\n\n";
        # print GATK "java $gtkMem -jar $GATK -R $hg19Fasta -T ReduceReads -log $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.table.printreads.datacompression.log -I $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.bam -o $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.reduced.bam\n\n";
	close GATK;
   }
}

system("mkdir -p $outputDir\/Variant");

###############################################################################
# Call variants with GATK UnifiedGenotyper
#
my $command2="java $gtkMem -jar $GATK -R $hg19Fasta -T UnifiedGenotyper -stand_call_conf $unifiedGenotyperStandcallconf -stand_emit_conf $unifiedGenotyperStandemitconf -glm $unifiedGenotyperGlm -out_mode $unifiedGenotyperOutmode";
for my $key (keys %SampleBamMap){
    $command2 .= " -I $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.bam";
}
open(VARIANT_UG, ">$outputDir\/variant.UnifiedGenotyper.sh");
print VARIANT_UG "\#!\/bin\/sh\n\n\n";
if ($parallel==1){
    open(Cluster, "<$cluster");
    while (<Cluster>){
	print VARIANT_UG $_;
    }
    close Cluster;
}
print VARIANT_UG "$command2 --dbsnp /media/disk2/hg19/dbsnp_138.hg19.vcf -o $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.raw.vcf >& $outputDir\/Variant\/UnifiedGenotyper.GATK.vcf.log\n\n";
if($VQSR==1){
    print VARIANT_UG "java $gtkMem -jar $GATK -T VariantRecalibrator -R $hg19Fasta -input $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.raw.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap -resource:omni,known=false,training=true,truth=false,prior=12.0 $omni -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp138 -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche [100.0, 99.9, 99.0, 90.0] -percentBad 0.01 -minNumBad 1000 -recalFile $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.raw.vcf.SNP.recal -tranchesFile $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.raw.vcf.SNP.tranches -rscriptFile $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.raw.vcf.SNP.plots.R\n\n";

    print VARIANT_UG "java $gtkMem -jar $GATK -T ApplyRecalibration -R $hg19Fasta -input $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.raw.vcf -mode SNP --ts_filter level 99.0 -recalFile $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.raw.vcf.SNP.recal -tranchesFile $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.raw.vcf.SNP.tranches -o $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.var.flt.vcf\n\n";
}else{
    print VARIANT_UG "java $gtkMem -jar $GATK -R $hg19Fasta -T VariantFiltration -o $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.var.flt.vcf --variant $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.raw.vcf --clusterSize 3 --clusterWindowSize 10 --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\" --filterExpression \"DP < 6\" --filterName \"LowCoverage\" --filterExpression \"QUAL < 30.0\" --filterName \"VeryLowQual\" --filterExpression \"QUAL >= 30.0 && QUAL < 50.0\" --filterName \"LowQual\" --filterExpression \"QD < 1.5\" --filterName \"LowQD\" --filterExpression \"SB > -10.0\" --filterName \"StrandBias\" >& $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.VariantFiltration.log\n\n";
}
print VARIANT_UG "java $gtkMem -jar $GATK -R $hg19Fasta -T SelectVariants -o $outputDir\/Variant\/allsamples.unifiedgenotyper.var.flt.excludeFiltered.vcf --variant $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.var.flt.vcf --excludeFiltered >& $outputDir\/Variant/allsamples.unifiedgenotyper.var.flt.vcf.excludeFiltered.log\n\n";
print VARIANT_UG "$convert2annovar -format vcf4old --includeinfo $outputDir\/Variant\/allsamples.unifiedgenotyper.var.flt.excludeFiltered.vcf > $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.var.flt.excludeFiltered.vcf.annovar 2>$outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.var.excludeFiltered.flt.vcf.annovar.log\n\n";
print VARIANT_UG "$summarizeannovar --buildver $buildver --verdbsnp $verdbsnp --ver1000g $ver1000g --veresp $veresp $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.var.flt.excludeFiltered.vcf.annovar $annovardatabase --outfile $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.var.flt.excludeFiltered.vcf.annovar.annotation >& $outputDir\/Variant\/allsamples.GATK.unifiedgenotyper.var.flt.excludeFiltered.vcf.annovar.annotation.log\n\n";
close VARIANT_UG;

###############################################################################
# Call variants with GATK HaplotypeCaller: 
#
my $command3="java $gtkMem -jar $GATK -R $hg19Fasta -T HaplotypeCaller -stand_call_conf $unifiedGenotyperStandcallconf -stand_emit_conf $unifiedGenotyperStandemitconf -minPruning 10";
for my $key (keys %SampleBamMap){
    $command3 .= " -I $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.bam";
}
open(VARIANT_HC, ">$outputDir\/variant.HaplotypeCaller.sh");
print VARIANT_HC "\#!\/bin\/sh\n\n\n";
if ($parallel==1){
    open(Cluster, "<$cluster");
    while (<Cluster>){
	print VARIANT_HC $_;
    }
    close Cluster;
}
print VARIANT_HC "$command3 --dbsnp /media/disk2/hg19/dbsnp_138.hg19.vcf -o $outputDir\/Variant\/allsamples.GATK.haplotypecaller.raw.vcf >& $outputDir\/Variant\/HaplotypeCallder.GATK.vcf.log\n\n";
if($VQSR==1){
    print VARIANT_HC "java $gtkMem -jar $GATK -T VariantRecalibrator -R $hg19Fasta -input $outputDir\/Variant\/allsamples.GATK.haplotypecaller.raw.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap -resource:omni,known=false,training=true,truth=false,prior=12.0 $omni -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp138 -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche [100.0, 99.9, 99.0, 90.0] -percentBad 0.01 -minNumBad 1000 -recalFile $outputDir\/Variant\/allsamples.GATK.haplotypecaller.raw.vcf.SNP.recal -tranchesFile $outputDir\/Variant\/allsamples.GATK.haplotypecaller.raw.vcf.SNP.tranches -rscriptFile $outputDir\/Variant\/allsamples.GATK.haplotypecaller.raw.vcf.SNP.plots.R\n\n";

    print VARIANT_HC "java $gtkMem -jar $GATK -T ApplyRecalibration -R $hg19Fasta -input $outputDir\/Variant\/allsamples.GATK.haplotypecaller.raw.vcf -mode SNP --ts_filter level 99.0 -recalFile $outputDir\/Variant\/allsamples.GATK.haplotypecaller.raw.vcf.SNP.recal -tranchesFile $outputDir\/Variant\/allsamples.GATK.haplotypecaller.raw.vcf.SNP.tranches -o $outputDir\/Variant\/allsamples.GATK.haplotypecaller.var.flt.vcf\n\n";
}else{
    print VARIANT_HC "java $gtkMem -jar $GATK -R $hg19Fasta -T VariantFiltration -o $outputDir\/Variant\/allsamples.GATK.haplotypecaller.var.flt.vcf --variant $outputDir\/Variant\/allsamples.GATK.haplotypecaller.raw.vcf --clusterSize 3 --clusterWindowSize 10 --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\" --filterExpression \"DP < 6\" --filterName \"LowCoverage\" --filterExpression \"QUAL < 30.0\" --filterName \"VeryLowQual\" --filterExpression \"QUAL >= 30.0 && QUAL < 50.0\" --filterName \"LowQual\" --filterExpression \"QD < 1.5\" --filterName \"LowQD\" --filterExpression \"SB > -10.0\" --filterName \"StrandBias\" >& $outputDir\/Variant\/allsamples.GATK.haplotypecaller.VariantFiltration.log\n\n";
}
print VARIANT_HC "java $gtkMem -jar $GATK -R $hg19Fasta -T SelectVariants -o $outputDir\/Variant\/allsamples.haplotypecaller.var.flt.excludeFiltered.vcf --variant $outputDir\/Variant\/allsamples.GATK.haplotypecaller.var.flt.vcf --excludeFiltered >& $outputDir\/Variant/allsamples.haplotypecallder.var.flt.vcf.excludeFiltered.log\n\n";
print VARIANT_HC "$convert2annovar -format vcf4old --includeinfo $outputDir\/Variant\/allsamples.haplotypecaller.var.flt.excludeFiltered.vcf > $outputDir\/Variant\/allsamples.haplotypecaller.var.flt.excludeFiltered.vcf.annovar 2>$outputDir\/Variant\/allsamples.GATK.haplotypecaller.var.flt.excludeFiltered.vcf.annovar.log\n\n";
print VARIANT_HC "$summarizeannovar --buildver $buildver --verdbsnp $verdbsnp --ver1000g $ver1000g --veresp $veresp $outputDir\/Variant\/allsamples.haplotypecaller.var.flt.excludeFiltered.vcf.annovar $annovardatabase --outfile $outputDir\/Variant\/allsamples.GATK.haplotypecaller.var.flt.excludeFiltered.vcf.annovar.annotation >& $outputDir\/Variant\/allsamples.GATK.haplotypecaller.var.flt.excludeFiltered.vcf.annovar.annotation.log\n\n";
close VARIANT_HC;

###############################################################################
# Call variants with samtools
#
my $command4="$samtools mpileup -SDug -f $hg19Fasta";
my $command5=" | /media/disk2/software/bin/bcftools view -bvcg - > $outputDir\/Variant\/allsamples.samtools.var.raw.bcf 2>$outputDir\/Variant\/allsamples.samtools.var.raw.bcf.log";
for my $key (keys %SampleBamMap){
    $command4 .= " $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.bam";
}
$command4 .=$command5;
open(VARIANT_ST, ">$outputDir\/variant.samtools.sh");
print VARIANT_ST "\#!\/bin\/sh\n\n\n";
if ($parallel==1){
    open(Cluster, "<$cluster");
    while (<Cluster>){
	print VARIANT_ST $_;
    }
    close Cluster;
}
print VARIANT_ST "$command4\n\n";
print VARIANT_ST "/media/disk2/software/bin/bcftools view $outputDir\/Variant\/allsamples.samtools.var.raw.bcf | vcfutils.pl varFilter -D1000 > $outputDir\/Variant\/allsamples.samtools.var.raw.vcf 2>$outputDir\/Variant\/allsamples.samtools.var.flt.vcf.log\n\n";
if($VQSR==1){
    print VARIANT_ST "java $gtkMem -jar $GATK -T VariantRecalibrator -R $hg19Fasta -input $outputDir\/Variant\/allsamples.samtools.raw.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap -resource:omni,known=false,training=true,truth=false,prior=12.0 $omni -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp138 -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche [100.0, 99.9, 99.0, 90.0] -percentBad 0.01 -minNumBad 1000 -recalFile $outputDir\/Variant\/allsamples.samtools.raw.vcf.SNP.recal -tranchesFile $outputDir\/Variant\/allsamples.samtools.raw.vcf.SNP.tranches -rscriptFile $outputDir\/Variant\/allsamples.samtools.raw.vcf.SNP.plots.R\n\n";

    print VARIANT_ST "java $gtkMem -jar $GATK -T ApplyRecalibration -R $hg19Fasta -input $outputDir\/Variant\/allsamples.samtools.raw.vcf -mode SNP --ts_filter level 99.0 -recalFile $outputDir\/Variant\/allsamples.samtools.raw.vcf.SNP.recal -tranchesFile $outputDir\/Variant\/allsamples.samtools.raw.vcf.SNP.tranches -o $outputDir\/Variant\/allsamples.samtools.var.flt.vcf\n\n";
}else{
    print VARIANT_ST "java $gtkMem -jar $GATK -R $hg19Fasta -T VariantFiltration -o $outputDir\/Variant\/allsamples.samtools.var.flt.vcf --variant $outputDir\/Variant\/allsamples.samtools.var.raw.vcf --clusterSize 3 --clusterWindowSize 10 --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\" --filterExpression \"DP < 6\" --filterName \"LowCoverage\" --filterExpression \"QUAL < 30.0\" --filterName \"VeryLowQual\" --filterExpression \"QUAL >= 30.0 && QUAL < 50.0\" --filterName \"LowQual\" --filterExpression \"QD < 1.5\" --filterName \"LowQD\" --filterExpression \"SB > -10.0\" --filterName \"StrandBias\"\n\n";
}
print VARIANT_ST "java $gtkMem -jar $GATK -R $hg19Fasta -T SelectVariants -o $outputDir\/Variant\/allsamples.samtools.var.flt.excludeFiltered.vcf --variant $outputDir\/Variant\/allsamples.samtools.var.flt.vcf --excludeFiltered >& $outputDir\/Variant/allsamples.samtools.var.flt.vcf.excludeFiltered.log\n\n";
print VARIANT_ST "$convert2annovar -format vcf4old --includeinfo $outputDir\/Variant\/allsamples.samtools.var.flt.excludeFiltered.vcf  > $outputDir\/Variant\/allsamples.samtools.var.flt.excludeFiltered.vcf.annovar 2>$outputDir\/Variant\/allsamples.samtools.var.flt.excludeFiltered.vcf.annovar.log\n\n";
print VARIANT_ST "$summarizeannovar --buildver $buildver --verdbsnp $verdbsnp --ver1000g $ver1000g --veresp $veresp $outputDir\/Variant\/allsamples.samtools.var.flt.excludeFiltered.vcf.annovar $annovardatabase --outfile $outputDir\/Variant\/allsamples.samtools.var.flt.excludeFiltered.vcf.annovar.annotation >& $outputDir\/Variant\/allsamples.samtools.var.flt.excludeFiltered.vcf.annovar.annotation.log\n\n";
close VARIANT_ST;

###############################################################################
# Call variants with SNVer:
# 
system("mkdir -p $outputDir\/Variant\/BAM");
open(VARIANT_SNVer, ">$outputDir\/variant.SNVer.sh");
print VARIANT_SNVer "\#!\/bin\/sh\n\n\n";
if ($parallel==1){
    open(Cluster, "<$cluster");
    while (<Cluster>){
	print VARIANT_SNVer $_;
    }
    close Cluster;
}
for my $key (keys %SampleBamMap){
    print VARIANT_SNVer "cp $outputDir\/$key\/$key\_merged.sorted.markdup.realigned.recal.bam $outputDir\/Variant/BAM\n\n";
};

print VARIANT_SNVer "java -jar /media/disk2/software/SNVerPool.jar -i $outputDir\/Variant/BAM -o $outputDir\/Variant/allsamples.SNVer -n 10 -r $hg19Fasta >& $outputDir\/Variant/allsamples.SNVer.log\n\n";
print VARIANT_SNVer "$bgzip $outputDir\/Variant/allsamples.SNVer.filter.vcf\n\n";
print VARIANT_SNVer "$tabix -p vcf $outputDir\/Variant/allsamples.SNVer.filter.vcf.gz\n\n";
print VARIANT_SNVer "$bgzip $outputDir\/Variant/allsamples.SNVer.indel.filter.vcf\n\n";
print VARIANT_SNVer "$tabix -p vcf $outputDir\/Variant/allsamples.SNVer.indel.filter.vcf.gz\n\n";
print VARIANT_SNVer "/media/disk2/software/vcftools_0.1.11/perl/vcf-concat $outputDir\/Variant/allsamples.SNVer.filter.vcf.gz $outputDir\/Variant/allsamples.SNVer.indel.filter.vcf.gz | gzip -c > $outputDir\/Variant/allsamples.SNVer.concat.vcf.gz\n\n";
print VARIANT_SNVer "/media/disk2/software/vcftools_0.1.11/perl/vcf-sort $outputDir\/Variant/allsamples.SNVer.concat.vcf.gz > $outputDir\/Variant/allsamples.SNVer.concat.sorted.vcf\n\n";
print VARIANT_SNVer "$convert2annovar -format vcf4old --includeinfo $outputDir\/Variant/allsamples.SNVer.concat.sorted.vcf > $outputDir\/Variant\/allsamples.SNVer.concat.sorted.vcf.annovar 2>$outputDir\/Variant\/allsamples.SNVer.concat.sorted.vcf.annovar.log\n\n";
print VARIANT_SNVer "$summarizeannovar --buildver $buildver --verdbsnp $verdbsnp --ver1000g $ver1000g --veresp $veresp $outputDir\/Variant\/allsamples.SNVer.concat.sorted.vcf.annovar $annovardatabase --outfile $outputDir\/Variant\/allsamples.SNVer.concat.sorted.vcf.annovar.annotation >& $outputDir\/Variant\/allsamples.SNVer.concat.sorted.vcf.annovar.annotation.log\n\n";
close VARIANT_SNVer;

###############################################################################
# Output consensus callset shared by different callers
# 
open(VARIANT4, ">$outputDir\/variant.summary.sh");
print VARIANT4 "\#!\/bin\/sh\n\n\n";
if ($parallel==1){
    open(Cluster, "<$cluster");
    while (<Cluster>){
	print VARIANT4 $_;
    }
    close Cluster;
}
#open(Cluster, "<$parallel");
#while (<Cluster>){
#print VARIANT4 $_;
#     }
#close Cluster;
print VARIANT4 "$bgzip $outputDir\/Variant\/allsamples.samtools.var.flt.excludeFiltered.vcf\n\n";
print VARIANT4 "$tabix -p vcf $outputDir\/Variant\/allsamples.samtools.var.flt.excludeFiltered.vcf.gz\n\n";
print VARIANT4 "$bgzip $outputDir\/Variant\/allsamples.haplotypecaller.var.flt.excludeFiltered.vcf\n\n";
print VARIANT4 "$tabix -p vcf $outputDir\/Variant\/allsamples.haplotypecaller.var.flt.excludeFiltered.vcf.gz\n\n";
print VARIANT4 "$bgzip $outputDir\/Variant\/allsamples.unifiedgenotyper.var.flt.excludeFiltered.vcf\n\n";
print VARIANT4 "$tabix -p vcf $outputDir\/Variant\/allsamples.unifiedgenotyper.var.flt.excludeFiltered.vcf.gz\n\n";
print VARIANT4 "$bgzip $outputDir\/Variant/allsamples.SNVer.concat.sorted.vcf\n\n";
print VARIANT4 "$tabix -p vcf $outputDir\/Variant/allsamples.SNVer.concat.sorted.vcf.gz\n\n";
print VARIANT4 "/media/disk2/software/vcftools_0.1.11/perl/vcf-isec -n +4 -f $outputDir\/Variant\/allsamples.samtools.var.flt.excludeFiltered.vcf.gz $outputDir\/Variant\/allsamples.haplotypecaller.var.flt.excludeFiltered.vcf.gz $outputDir\/Variant\/allsamples.unifiedgenotyper.var.flt.excludeFiltered.vcf.gz $outputDir\/Variant/allsamples.SNVer.concat.sorted.vcf.gz | bgzip -c > $outputDir\/Variant\/consensus.vcf.gz\n\n";
# annovar
print VARIANT4 "$convert2annovar -format vcf4old --includeinfo $outputDir\/Variant\/consensus.vcf > $outputDir\/Variant\/consensus.vcf.annovar 2>$outputDir\/Variant\/consensus.vcf.annovar.log\n\n";
print VARIANT4 "$summarizeannovar --buildver $buildver --verdbsnp $verdbsnp --ver1000g $ver1000g --veresp $veresp $outputDir\/Variant\/consensus.vcf.annovar $annovardatabase --outfile $outputDir\/Variant\/consensus.vcf.annovar.annotation >& $outputDir\/Variant\/consensus.vcf.annovar.annotation.log\n\n";
# VEP Ensembl
print VARIANT4 "perl $VEP -i $outputDir\/Variant\/consensus.vcf -o $outputDir\/Variant\/consensus.vcf.VEP --offline --cache";
close VARIANT4;

