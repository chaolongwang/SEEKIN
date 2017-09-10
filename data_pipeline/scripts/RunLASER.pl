
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
my @DEPENDS=("Foundation")        ;
my $VERSION = "Version 1.00"        ;
my $default_str="this is a default string";
my $MAX = 10**(9);
my $app = "$Bin/../apps";
my %global=();
my %centpromer=();
my $beagleMerge="";
sub setDefault();
sub setParametersFromCommandLine();
sub setParametersFromFile($@);
sub getPileup2SeqCMD();


setDefault();
setParametersFromCommandLine();
getPileup2SeqCMD();




sub getPileup2SeqCMD(){

  my $BATCH_SIZE=$global{BATCH_SIZE};
  my $bamList=$global{BAM_LIST};
  my $samtoolsPath=$global{SAMTOOLS};
  my $refPath=$global{REF_FILE};
  my $refBed=$global{REFBED_FILE};
  my $refSite=$global{REFSITE_FILE};
  my $geno=$global{REFGENO_FILE};
  my $LASER=$global{LASER};
  my $qsub=$global{QSUB};
  my $coord=$global{REFCOORD_FILE};
  my $pileup2seq=$global{PILEUP2SEQ};
  my $OUTPUT_DIR=abs_path($global{OUTPUT_DIR});
  my $output=$OUTPUT_DIR."/".$global{OUT_DIR};
  my $SEEKIN=$global{SEEKIN};
  my $BETA_FILE=$global{BETA_FILE};
  my $libraryID=$global{SWAP_ID};
  my $coordList="";
  open OO, ">$OUTPUT_DIR/jobfiles/runLASER_batch.sh";
  open O1, ">$output/ID.index";

  mkdir($output, 0700) unless(-d $output);
  my $N = 0;
  open (FIN, $bamList) || die "can't open $global{BAM_LIST}: $!";
  while(my $line=<FIN>){ $N++; }
  close FIN;
  print "$N bam files in total.\n";


  
  my $id="";
  $BATCH_SIZE=39;
  open (FIN, $bamList) || die "can't open $bamList: $!";
  open T, ">$OUTPUT_DIR/tmp/laser.batch.lst";
  for(my $j=0; $j<($N/$BATCH_SIZE); $j++){
          print T "$j\n";
          my $batchfile = "$OUTPUT_DIR/jobfiles/RunLASER$j.sh";
          my $str="python $pileup2seq \\\n  -m $refSite \\\n  -o $output/pileup2seq$j";
          open (OUT, ">$batchfile") || die "can't open $batchfile: $!";
          for(my $i=0; $i<$BATCH_SIZE; $i++){
                  my $bam = <FIN>;
                  chomp($bam);    
                  my @f=split(/\//,$bam);   
                  $id = $f[scalar(@f)-1];
                  my @tmp=split(/\./,$id);
                  print O1 "$j"."_$i $tmp[0]\n";
                  print OUT "$samtoolsPath mpileup -q 20 -Q 20 \\\n  -f $refPath  \\\n  -l $refBed \\\n  $bam \\\n  >$output/$j"."_$i.pileup\n\n";
          $str=$str."\\\n  $output/$j"."_$i.pileup ";
                  
                  last if eof(FIN);
          }
          $str=$str."\\\n";
          print OUT "$str\n\n";
          print OUT "$global{LASER}  \\\n  -g $geno   \\\n  -s $output/pileup2seq$j.seq   \\\n  -c $coord \\\n  $global{LASER_PAR_SET} \\\n  -o $output/$j\n\n";
          close OUT;
          print OO "$qsub runLASER$j -pe OpenMP 1 bash $batchfile\n";
  }
  close FIN;
  close O1;
  close T;
  my @t=split(/\s+|\t/,$global{LASER_PAR_SET});
  my $pca_seekin=2;

  for(my $l=0;$l<@t;$l++){
	   if($t[$l] eq "-k"){$pca_seekin=$t[$l+1];}
  }
  open O3, ">$OUTPUT_DIR/jobfiles/runSEEKIN.sh\n";
  open O2, ">$OUTPUT_DIR/jobfiles/mergeLaserOutput.sh\n";
  print O2 "cat laser/*.SeqPC.coord  |  awk '!/popID/ || NR==1'  > $output/laser.seqPC.coord\n";
  print O3 "\n# Generate the beta file for the chosen refernece panel\n\n";
  print O3 "$SEEKIN modelAF  \\\n  -i $global{LASER_REF_VCF}  \\\n  -c $coord  \\\n  -k  $pca_seekin  \\\n  -o $output/REF.beta\n\n";
  print O3 "$SEEKIN  getAF  \\\n  -i  $output/laser.seqPC.coord -k $pca_seekin \\\n  -b  $output/REF.beta  \\\n  -o  $output/All.af\n\n";
  print O3 "bcftools reheader -s $output/ID.index   $output/All.af.gz  \\\n  > $output/seekin.final.af.gz \n\n";
  print O3 "tabix -p vcf $output/seekin.final.af.gz \n\n";
  print O3 "$SEEKIN  kinship  \\\n  $global{SEEKIN_PAR_SET} \\\n  -i ./snp/Beagle.gp.vcf.gz \\\n  -f ./laser/seekin.final.af.gz\\\n  -o ./seekin/seekin \n";
  close O2;
  close O3;
}


sub setParametersFromCommandLine(){

    my @specOpts =() ;
    my $result = $tf->TIGR_GetOptions(
       "c=s"                   => \$global{"CONFIGURE_FILE"},
       "i=s"                   => \$global{"BAM_LIST"},
       "r=s"                   => \$global{"REF_FILE"},
       "g=s"                   => \$global{"REFBED_FILE"},
       "s=s"                   => \$global{"REFSITE_FILE"},
       "p=s"                   => \$global{"PILEUP2SEQ"},
       "d=s"                   => \$global{"REFCOORD_FILE"},
       "l=i"                   => \$global{"LASER"},
       "m=i"                   => \$global{"SAMTOOLS"},
       "t=i"                   => \$global{"LASER_BATCH_SIZE"},
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

    #print STDERR "\n";
    #print STDERR "###\n";
    #print STDERR "###  Read options from '$specFile' and generate the job files\n";
    #print STDERR "###  Please wait for a few minutes ... \n";
    #print STDERR "\n";

    open(F, "< $specFile") or die("Couldn't open '$specFile'", undef);

    while (<F>) {
        s/^\s+//;
        s/\s+$//;
        next if (m/^\s*\#/);
        next if (m/^\s*$/) ;
        if($_=~/option/){
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
    $global{REFGENO_FILE}="$global{LASER_REF_PREFIX}".".geno";
    $global{REFCOORD_FILE}="$global{LASER_REF_PREFIX}".".RefPC.coord";
    $global{REFBED_FILE}="$global{LASER_REF_PREFIX}".".bed";
    $global{REFSITE_FILE}="$global{LASER_REF_PREFIX}".".site";
    $global{SAMTOOLS}="$app/samtools";
    $global{SEEKIN}="$app/../../bin/seekin";
    $global{LASER}="$app/LASER/laser";
    $global{PILEUP2SEQ}="$app/LASER/pileup2seq.py";
}



sub setDefault(){

  my $HELPTEXT = qq~
  Program: runLASER (run LASER from WES datasets)
  Version: 1.000

  Usage:  perl RunLASER.pl [options] 

  Options: -c FILE  configure file 
     -i FILE  list of input BAM filenames, one per line [null]
     -b FILE  list of positions (chr pos) or regions (BED) to remove [null]
     -t INT   number of threads [1]
     -o STR   output VCF file prefix
     -h     display the help information
  ~;


my $MOREHELP = qq~
Return Codes:   0 - on success, 1 - on failure.
~;


    $global{"HELP"}=$HELPTEXT; 
    $global{"CONFIGURE_FILE"}=undef,
    $global{"BAM_LIST"}=undef,
    $global{"REF_FILE"}=undef,
    $global{"REFBED_FILE"}=undef,
    $global{"REFSITE_FILE"}=undef,
    $global{"PILEUP2SEQ"}=undef,
    $global{"REFCOORD_FILE"}=undef,
    $global{"LASER"}=undef,
    $global{"SAMTOOLS"}=undef,
    $global{"LASER_BATCH_SIZE"}=undef,
    $global{"SEEKIN_REF_VCF"} = undef,
    $global{"OUT_DIR"}="laser",
    $global{"QSUB"}="qsub -b y  -l mem_free=10G,h_rt=72:0:0  -e ./logs   -o  ./logs   -q long.q  -v PATH -cwd -terse -N";


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
