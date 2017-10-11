#!/usr/bin/perl
###############################################################################
#   Copyright (C) 2016, Genome Institute of Singapore (GIS);
#   all rights reserved. Authored by: Jinzhuang Dou
#   douj@gis.a-star.edu.sg
#
#
#   Redistribution and use in source and binary forms, with 
#   or without modification, are permitted provided that the 
#   following conditions are met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the 
#     documentation and/or other materials provided with the distribution
#
#   * Neither the name of the OUC nor the names of its contributors may be 
#     used to endorse or promote products derived from this software 
#     without specific prior written permission.
#
#############################################################################

use strict;
use warnings;
use Cwd 'abs_path';
use Time::localtime;
use FindBin qw($Bin); 
use lib "$Bin/../lib";
use Foundation;
my $tf = new Foundation                 ;
my $PRG = $tf->getProgramInfo('name')   ;
my $REV="1.0"                           ;                   
my @DEPENDS=("Foundation")				;
my $VERSION = "Version 1.00"  			;
my $default_str="this is a default string";
my $MAX = 10**(9);
my $app = "$Bin/../apps";
my $resource = "$Bin/../resource";
my %global=();
my %centpromer=();
my $beagleMerge="";
sub setDefault();
sub writeCMD;
sub getCentromerePos;
sub setParametersFromCommandLine();
sub setParametersFromFile($@);
sub getCMD();
sub readFileList($);


setDefault();
setParametersFromCommandLine();


if(defined $global{"CENTROMERE_FILE"}  && -e $global{"CENTROMERE_FILE"} ) {
	%centpromer = getCentromerePos( $global{"CENTROMERE_FILE"}, 1);
}
else{
	print "\nWarnings: no centromere position file is specified, we will ignore it.\n\n";
	for(my $i=1;$i<=22;$i++){
		$centpromer{$i}{"start"}=$MAX;
		$centpromer{$i}{"end"}=$MAX;
	}
}
$global{OUTPUT_DIR}=abs_path($global{OUTPUT_DIR});

getCMD();


sub getCMD(){

	my $folder="$global{OUTPUT_DIR}/jobfiles";
	my $logDir="$global{OUTPUT_DIR}/log";
	my $tmpDir ="$global{OUTPUT_DIR}/tmp";
	my %chr_freq_list=();
	my %chr_vcf_list=();
	my %chr_bed_list=();

	getSiteInfo($global{VCF_SITE_FILE});

	$global{REF_FREQ_LIST}="$tmpDir/sites.bed.lst";
	$global{REF_BED_LIST}="$tmpDir/sites.bed.lst";
	
	
	
	%chr_freq_list=readFileList($global{REF_FREQ_LIST});
	%chr_vcf_list=readFileList($global{BEAGLE_REF_LIST});
	%chr_bed_list=readFileList($global{REF_BED_LIST});

	my $chr_freq_cnt=scalar keys %chr_freq_list;
	my $chr_vcf_cnt=scalar keys %chr_vcf_list;
	my $chr_bed_cnt=scalar keys %chr_bed_list;



	if($chr_freq_cnt!=$chr_vcf_cnt || $chr_vcf_cnt!=$chr_bed_cnt){
		print "The CHR number in $global{REF_FREQ_LIST} is $chr_freq_cnt\n";
		print "The CHR number in $global{BEAGLE_IMPUTE_FILE_LST} is $chr_vcf_cnt\n";
		print "The CHR number in $global{REF_BED_LIST} is $chr_bed_cnt\n";
		print "Please check them!\n";
		exit(-1);
	}
	else{
		my $ii=0;
		my $outputID = "";
		my $chuck_size_default = $global{"CHUNCK_SIZE"};
		open O1, ">$tmpDir/region.Lst";
		for($ii=1;$ii<=$chr_bed_cnt;$ii++){
			my $line1=0;
			my $line2=0;
			my %hash1=();
			my %hash2=();
			my %chrID=();
			my $hold_id="";

			$beagleMerge="";

			print "Generate the command [chr$ii]...\n";
			open F, "$chr_bed_list{$ii}";
			while(<F>){
					chomp;
					my @f=split(/\s+|\t/,$_);
					if($f[0]=~/chr/){$f[0]=substr($f[0],3);}
					if($f[1]<$centpromer{$f[0]}{"start"}){
						$line1++;
						$hash1{$line1}=$f[1];
						$chrID{$line1}=$f[0];
					}
					elsif ($f[1]>$centpromer{$f[0]}{"end"}){
						$line2++;
						$hash2{$line2}=$f[1];
						$chrID{$line2}=$f[0];
					}
			}
			close F;
			####
			####  modify the CHUNCK_SIZE when the marker number is smaller than specified chunck size
			####  


			mkdir($global{OUT_DIR},0700) unless (-d $global{OUT_DIR} );
			mkdir($folder, 0700) unless(-d $folder);
			
			if($line1<$chuck_size_default) {$global{CHUNCK_SIZE} = $line1;}
			else{ $global{CHUNCK_SIZE} = $chuck_size_default;}
	
			for(my $i=1;$i<=$line1-$global{"CHUNCK_SIZE"}/2;$i=$i+$global{"CHUNCK_SIZE"}){
				my $tmp=$i+$global{"CHUNCK_SIZE"}+$global{"CHUNCK_OVERLAP"};
				if ($line1-$tmp<=$global{"CHUNCK_SIZE"}/2){$tmp=$line1;}
				open O, ">$folder/getGP_chr"."$chrID{$i}"."_"."$hash1{$i}"."_"."$hash1{$tmp}".".sh";
				my $region="$chrID{$i}:$hash1{$i}-$hash1{$tmp}";
				my $outputID="$chrID{$i}_$hash1{$i}_$hash1{$tmp}";
				my $cmd=writeCMD($region, $outputID, $chr_bed_list{$ii}, $chr_freq_list{$ii}, $chr_vcf_list{$ii});
				print O "$cmd\n";
				close O;
				print O1 "$chrID{$i}"."_"."$hash1{$i}"."_"."$hash1{$tmp}\n";
			}
			
			if($line2<$chuck_size_default) {$global{CHUNCK_SIZE} = $line2;}
			else { $global{CHUNCK_SIZE} = $chuck_size_default;}

			for(my $i=1;$i<=$line2-$global{"CHUNCK_SIZE"}/2;$i=$i+$global{"CHUNCK_SIZE"}){
				my $tmp=$i+$global{"CHUNCK_SIZE"}+$global{"CHUNCK_OVERLAP"};
				if ($line2-$tmp<=$global{"CHUNCK_SIZE"}/2){$tmp=$line2;}
				open O, ">$folder/getGP_chr"."$chrID{$i}"."_"."$hash2{$i}"."_"."$hash2{$tmp}".".sh";
				my $region="$chrID{$i}:$hash2{$i}-$hash2{$tmp}";
				$outputID="$chrID{$i}_$hash2{$i}_$hash2{$tmp}";  # Use this variable latter
				my $cmd=writeCMD($region, $outputID, $chr_bed_list{$ii}, $chr_freq_list{$ii}, $chr_vcf_list{$ii});
				print O "$cmd\n";
				close O;
				print O1 "$chrID{$i}"."_"."$hash2{$i}"."_"."$hash2{$tmp}\n";
			}
			open O, ">$folder/merge_chr$ii".".sh";
			$beagleMerge=~s/^ +//; 
			my @tmp = split(/\s+|\t/,$beagleMerge);
			if (scalar(@tmp) > 1) {
				my $mergeCMD="$global{JAVA}  -jar $app/BEAGLE/mergevcf.jar  $ii  $beagleMerge > $global{OUT_DIR}/chr$ii.gp.vcf";
				print O "$mergeCMD\n";
				$mergeCMD="bgzip  -f  $global{OUT_DIR}/chr$ii.gp.vcf";
				print O "$mergeCMD\n";
			}
			else{
				my $mergeCMD="cp $global{OUT_DIR}/$outputID.gp.vcf.gz  $global{OUT_DIR}/chr$ii.gp.vcf.gz";
				print O "$mergeCMD\n";
				
			}
			close O;
			#print OO "qsub -b y  -l mem_free=20G,h_rt=2:0:0  -e ./log   -o  ./log   -q short.q  -v PATH -cwd -terse -N snpChrMerge$ii -hold_jid $hold_id  -pe OpenMP 1 bash $folder/merge_chr$ii".".sh\n\n";
			#print O1 "qsub -b y  -l mem_free=20G,h_rt=2:0:0  -e ./log   -o  ./log   -q short.q  -v PATH -cwd -terse -N snpChrMerge$ii -pe OpenMP 1 bash $folder/merge_chr$ii".".sh\n\n";
		}
		close O1;
	}

	open O, ">$folder/merge_chr.sh"; 
	my $merge_hold_id="snpChrMerge1";
	my $merge_chr_str="$global{OUT_DIR}/chr1.gp.vcf.gz";
	for(my $i=2;$i<=22;$i++){
		$merge_hold_id=$merge_hold_id.",snpChrMerge$i";
		$merge_chr_str=$merge_chr_str." $global{OUT_DIR}/chr$i.gp.vcf.gz";
	}
	## merge vcf files from different chr other than using vcf-concat
	print O "zcat $merge_chr_str | sed '500,1000000000000{/#/d;}'  |  bgzip -c > $global{OUT_DIR}/Beagle.gp.vcf.gz\n\n";
	close O;
	#print OO "qsub -b y  -l mem_free=10G,h_rt=2:0:0  -e ./log   -o  ./log   -q short.q  -v PATH -cwd -terse -N MergeChr -hold_jid $merge_hold_id  -pe OpenMP 1 bash $folder/merge_chr.sh\n\n"

}


sub writeCMD{
	my $CMD="step=\$1\n";
	my $region=shift;
	my $outputID=shift;
	my $bedfile=shift;
	my $freqfile=shift;
	my $vcffile=shift;
	$CMD=$CMD." "."if [ \$step -le 1 ];then\n$global{SAMTOOLS} mpileup -b $global{BAM_LIST} -l $bedfile -f $global{REF_FILE} -r $region -q 20 -Q 20 -t DP -v  -o $global{OUT_DIR}/$outputID.gl.vcf.gz\nfi\n";
	$CMD=$CMD." "."if [ \$step -le 2 ];then\n$global{VCFORMAT}  -i  $global{OUT_DIR}/$outputID.gl.vcf.gz  -o $global{OUT_DIR}/$outputID.format.vcf.gz -f $freqfile\nfi\n";
        
	if($global{BEAGLE_IMPUTE_MODE}==1){
		$CMD=$CMD."if [ \$step -le 3 ];then\n$global{JAVA}  -jar $global{BEAGLE}  gl=$global{OUT_DIR}/$outputID.format.vcf.gz  ref=$vcffile  chrom=$region  out=$global{OUT_DIR}/$outputID.gp  $global{BEAGLE_PAR_SET}\nfi\n";
    }
    elsif($global{BEAGLE_IMPUTE_MODE}==0){
    	$CMD=$CMD."if [ \$step -le 3 ];then\n$global{JAVA}  -jar $global{BEAGLE}  gl=$global{OUT_DIR}/$outputID.format.vcf.gz  chrom=$region  out=$global{OUT_DIR}/$outputID.gp  $global{BEAGLE_PAR_SET}\nfi\n";

    }
    $beagleMerge=$beagleMerge."  $global{OUT_DIR}/$outputID.gp.vcf.gz";
    return $CMD;
}


sub readFileList($){

	my $inputfile =shift @_;
	my %chr_list=();
	my $chrCNT=0;
	if(!-e $inputfile) {
		print "Error: no $inputfile is specified, please check it!\n";
		$tf  -> printUsageInfoAndExit();
	}
	else{
		open FIN, "<$inputfile";
		while(<FIN>){
			chomp;
			if(length($_)>1){
					$chrCNT++;
					$chr_list{$chrCNT}=$_;
			}
		}
		close FIN;
	}
	return %chr_list;
}

sub getSiteInfo{
	my $inputfile = shift;
	if(!-e $inputfile) {
		print "Error: no $inputfile is specified, please check it!\n";
		$tf  -> printUsageInfoAndExit();
	}
	else{
		if ($inputfile =~ m/\.gz$/){open(IN,"gunzip -c $inputfile |");}
		else{open(IN,$inputfile) || die "Cannot open $inputfile\n";}
		my $tmpDir = "$global{OUTPUT_DIR}/tmp";
		system("mkdir -p $tmpDir");
		my $chr_pre = "Null";
		my $str = "";
		
		while(<IN>){
			chomp;
			# if($_=~/#/ || $_=~/overlap/){;}  # remove the overlap flag 
			if($_=~/#/){;}  
			else{
				my @f=split(/\s+|\t/,$_);
				if ($f[0] eq "$chr_pre"){
					my $start = $f[1] -1;
					print O "$f[0]\t$start\t$f[1]\t$f[3]\t$f[4]\n";
				}
				else{
					close(O);
					$str = "$str"."$tmpDir/chr$f[0].sites.bed\n";
					open O, ">$tmpDir/chr$f[0].sites.bed";
					my $start = $f[1] -1;
					print O "$f[0]\t$start\t$f[1]\t$f[3]\t$f[4]\n";
					$chr_pre = $f[0];
				}
		    }
		}
		close(IN);
		close O;
		open O, ">$tmpDir/sites.bed.lst";
		print O "$str";
		close O;
	}

}
sub getCentromerePos{
	my $inputfile =shift;
	my %hash=();

	if (!-e $inputfile) { print "ERROR, the $inputfile is not identified, please check it!\n";}
	else{
		open F, "<$inputfile";
		while(<F>){
			chomp;
			if($_=~/acen/)
			{
				my @f=split(/\s+|\t/,$_);
				if($f[0]=~/chr/){$f[0]=substr($f[0],3);}
				if(exists $hash{$f[0]}){
					$hash{$f[0]}{"end"}=$f[2];
				}
				else{
					$hash{$f[0]}{"start"}=$f[1];
				}
			}
		}
		close F;
		return %hash;
	}
}



sub setParametersFromCommandLine(){

    my @specOpts =() ;
    my $result = $tf->TIGR_GetOptions(
	   "c=s"                   => \$global{"CONFIGURE_FILE"},
       "i=s"                   => \$global{"BAM_LIST"},
       "N=i"                   => \$global{"CHUNCK_SIZE"},
       "n=i"                   => \$global{"CHUNCK_OVERLAP"},
       "o=s"				   => \$global{"OUTPUT_DIR"}
    );

    if (@specOpts || !$result) {
       $tf  -> printUsageInfoAndExit();
    }
    if(defined $global{"CONFIGURE_FILE"}) {
       setParametersFromFile($global{"CONFIGURE_FILE"});   
    }
}

sub setParametersFromFile($@) {


    my $specFile  = shift @_;

    #  Client should be ensuring that the file exists before calling this function.
    die "Error: CONFIGURE_FILE '$specFile' not found.\n"  if (! -e "$specFile");

    print STDERR "\n";
    print STDERR "###\n";
    print STDERR "###  Read options from '$specFile' and generate the job files\n";
    print STDERR "###  Please wait for a few minutes ...";
    print STDERR "\n";

    open(F, "< $specFile") or die("Couldn't open '$specFile'", undef);

    while (<F>) {
        s/^\s+//;
        s/\s+$//;
        next if (m/^\s*\#/);
        next if (m/^\s*$/) ;
	if($_=~/JAVA/){
		my ($var, $val,$par) = split /\s+/;
		$val = "$val $par";
		setGlobal($var, $val);
	}
        elsif($_=~/option/){
          my ($var, $val, $tmp) = split /"/;
          my @t = split (/\s+/,$var);
          setGlobal($t[0], $val);
        }
	else{
        	my ($var, $val) = split /\s+/;
        	setGlobal($var, $val);
	}
    }
    close(F);
}



sub setDefault(){

	my $HELPTEXT = qq~
	Program: getGP (Tools for getting the genotype probabilities from low-coverage seqeuncing reads)
	Version: 1.000

	Usage:	perl getGP.pl [options] 

	Options: -c FILE  configure file 
		 -i FILE  list of input BAM filenames, one per line [null]
		 -N INT   number of VCF records per region [10,000]
		 -n INT   number of VCF records shared between consecutive regions [1,000]
		 -t INT   number of threads [1]
		 -h 	  display the help information
	~;


my $MOREHELP = qq~
Return Codes:   0 - on success, 1 - on failure.
~;


    $global{"HELP"}               = $HELPTEXT; 
    $global{"BAM_LIST"}           = undef    ;
    $global{"VCF_SITE_FILE"}      = undef    ;
    $global{"CENTROMERE_FILE"}    = "$resource/cytoBand.txt"    ;
    $global{"CONFIGURE_FILE"}     = undef    ;
    $global{"CHUNCK_SIZE"}        = 10000    ;
    $global{"CHUNCK_OVERLAP"}     = 1000     ;
    $global{"OUT_DIR"}            = "snp"    ;
    $global{"REF_FILE"}           = undef    ;
    $global{"SAMTOOLS"}           = "$app/samtools"    ;
    $global{"BEAGLE"}        	  = "$app/BEAGLE/beagle.27Jul16.86a.jar"  ;
    $global{"VCFORMAT"}			  = "$app/vcfFormat";
    $global{"BEAGLE_THREAD_NUM"}  = 5  ;
    $global{"BEAGLE_IMPUTE_MODE"} = 1  ;
    $global{"SEEKIN"}        = undef    ;
    $global{"MODE"}               = "run"    ;



    $tf->setHelpInfo($HELPTEXT.$MOREHELP);
    $tf->setUsageInfo($HELPTEXT);
    $tf->setVersionInfo($REV);
    $tf->addDependInfo(@DEPENDS);


}


sub getGlobal ($) {
    my $var = shift @_;
    if (!exists($global{$var})) {
      return undef;
    }
    return($global{$var});
}

sub setGlobal ($) {
    my $var = shift @_;
    my $val = shift @_;
    # $val = undef  if ($val eq "");
    $global{$var} = $val;
}
