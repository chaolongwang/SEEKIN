use Switch;
use File::Basename;

$dir=dirname(__FILE__);
$plotPath=$dir;
$true=$ARGV[0];
$test=$ARGV[1];
$format=$ARGV[3];
$output=$ARGV[2];
$pop=$ARGV[4];

print "true   $true\n";
print "test   $test\n";
print "format $format\n";
print "output $output\n";

open F, "</mnt/projects/wangcl/merck/kinship/MerckTest/KinshipTestSet/admixture/Merck/result_28031016/sampleInfo/CH.list";
$popCHS={};
$popMAS={};
while(<F>){
	chomp;
	@f=split(/\s+|\t/,$_);
	$pop{$f[0]}=1;
}
close F;

open F, "</mnt/projects/wangcl/merck/kinship/MerckTest/KinshipTestSet/admixture/Merck/result_28031016/sampleInfo/MA.list";
while(<F>){
	chomp;
	@f=split(/\s+|\t/,$_);      
	$pop{$f[0]}=2;
}         
close F;	
						
							$line=0;
open F, "<$true";
while(<F>){
	chomp;
	$line++;
	@f=split(/\s+|\t/,$_);
	if($line>1){
		$kin_true{$f[0]}{$f[1]}="$f[2]\t$f[3]\t$f[4]";
		$kin_true{$f[1]}{$f[0]}="$f[2]\t$f[3]\t$f[4]";
	}
}
close F;
close F;

$line=0;
print "$format\n";
switch ($format) {
      case "SEEKIN"  {getFormatSEEKIN();}
      case "PCrelate" {getFormatPCrelate();}
      case "REAP" {getFormatREAP();}
      case "KING-robust" {getFormatKINGrobust();}
      case "KING-homo" {getFormatKINGhomo();}
      case "lcmLkin"  { print "Yes\n"; getFormatlcmLkin(); }
      case "GCTA"  {getFormatGCTA();}
      case "GRM" {getFormatGRM();}
      case "RelateAdmix" {getFormatRelateAdmix();}
}

if($pop eq "homo"){
      system("Rscript $plotPath/plot.r  $output.plot.fmt  $output");
}
else{
      system("Rscript $plotPath/plot.r  $output.plot.fmt  $output");
}

sub getRelateAdmix(){
	$input=$test;
}

sub getFormatSEEKIN(){
	$input=$test;
	open F, "<$input";
	open O, ">$output".".plot.fmt";
	while(<F>){
		chomp;
		$line++;
		@f=split(/\s+|\t/,$_);
		$POP=-1;
		if($line>1){
			if( exists $kin_true{$f[0]}{$f[1]} ){

				print O "$_\t0\t$kin_true{$f[0]}{$f[1]}\t$pop{$f[0]}$pop{$f[1]}\n";
	    	}
			#else{  print "Not find paris $f[0]<->$f[1] in $input_lcmlkin\n";}
		}
		else{
			print O "Ind1\tInd2\tNSNP\tKin_est\tIBD0_est\trelationship\tKin_true\tIBD_true\tclass\n";

		}
	}
	close F;
	close O;
}


sub getFormatPCrelate(){
	$input=$test;
	open F, "<$input";
	open O, ">$output".".plot.fmt";
	while(<F>){
		chomp;
		$line++;
		@f=split(/\s+|\t/,$_);
		$POP=-1;
		if($line>1){
			if( exists $kin_true{$f[2]}{$f[1]} ){

				print O "$f[0]\t$f[2]\t$f[4]\t$f[8]\t$f[5]\t$kin_true{$f[2]}{$f[1]}\t$pop{$f[2]}$pop{$f[1]}\n";
	    	}
			#else{  print "Not find paris $f[0]<->$f[1] in $input_lcmlkin\n";}
		}
		else{
			print O "Ind1\tInd2\tNSNP\tKin_est\tIBD0_est\trelationship\tKin_true\tIBD0_true\tclass\n";

		}
	}
	close F;
	close O;
}

sub getFormatlcmLkin(){
	$input=$test;
	print "Yes, you are using lcmLkin mode\n";
	open F, "<$input";
	open O, ">$output".".plot.fmt";
	while(<F>){
		chomp;
		$line++;
		@f=split(/\s+|\t/,$_);
		$POP=-1;
		if($line>1){
			if(exists $kin_true{$f[0]}{$f[1]} ){
				$kin_est=$f[5]/2;
				print O "$f[0]\t$f[1]\t$f[6]\t$kin_est\t$f[2]\t$kin_true{$f[0]}{$f[1]}\t$pop{$f[0]}$pop{$f[1]}\n";
			}
		}
		else{
				print O "Ind1\tInd2\tNSNP\tKin_est\tIBD0_est\trelationship\tKin_true\tIBD0_true\tclass\n";
 				#print O "Ind1\tInd2\tNSNP\tKin_est\tIBD0_est\trelationship\tKin_true\tIBD0_true\tclass\n";
		}
	}
	close F;
	close O;
}

sub getFormatGRM(){
	$input=$test;
 	open F, "<$input".".grm.id";
        open O, ">$output".".plot.fmt";

	$line=0;
        while(<F>){
 	  $line++;
          @f=split(/\s+|\t/,$_);
          $list{$line}=$f[0];
        }
        close F;
        print "$input".".grm.gz\n";
        open F, "<$input".".grm.gz";
        print O "Ind1\tInd2\tNSNP\tKin_est\tIBD0_est\trelationship\tKin_true\tIBD0_true\tclass\n";
        while(<F>){
           @f=split(/\s+|\t/,$_);
           $f[0]=$list{$f[0]};
           $f[1]=$list{$f[1]};
           $f[3]=$f[3]/2;
           if( exists $kin_true{$f[0]}{$f[1]} ){
               print O "$f[0]\t$f[1]\t$f[3]\t$f[3]\t$f[3]\t$kin_true{$f[0]}{$f[1]}\t$pop{$f[0]}$pop{$f[1]}\n";
           }
        } 
        close O;
        close F;
}


sub getFormatREAP(){
	$input=$test;
	open F, "<$input";
	open O, ">$output".".plot.fmt";
	while(<F>){
		chomp;
		$line++;
		@f=split(/\s+|\t/,$_);
		$POP=-1;
		if($line>1){
			if( exists $kin_true{$f[2]}{$f[1]} ){

				print O "$f[0]\t$f[2]\t$f[4]\t$f[8]\t$f[5]\t$kin_true{$f[2]}{$f[1]}\t$pop{$f[2]}$pop{$f[1]}\n";
	    	}
			#else{  print "Not find paris $f[0]<->$f[1] in $input_lcmlkin\n";}
		}
		else{
			print O "Ind1\tInd2\tNSNP\tKin_est\tIBD0_est\trelationship\tKin_true\tIBD0_true\tclass\n";

		}
	}
	close F;
	close O;
}


sub getFormatGCTA(){
	$input=$test;
	open F, "<$input";
	open O, ">$output".".plot.fmt";
	while(<F>){
		chomp;
		$line++;
		@f=split(/\s+|\t/,$_);
		$POP=-1;
		if($line>1){
			if( exists $kin_true{$f[0]}{$f[1]} ){
				$f[2]=$f[2]/2;
				print O "$f[0]\t$f[1]\t0\t$f[2]\t0\t$kin_true{$f[0]}{$f[1]}\t$pop{$f[0]}$pop{$f[1]}\n";
	    	}
			#else{  print "Not find paris $f[0]<->$f[1] in $input_lcmlkin\n";}
		}
		else{
			print O "Ind1\tInd2\tNSNP\tKin_est\tIBD0_est\trelationship\tKin_true\tIBD0_true\tclass\n";

		}
	}
	close F;
	close O;
}








sub getFormatKINGhomo(){
	$input=$test;
	open F, "<$input";
	open O, ">$output".".plot.fmt";
	while(<F>){
		chomp;
		$line++;
		@f=split(/\s+|\t/,$_);
		$POP=-1;
		if($line>1){
			if( exists $kin_true{$f[2]}{$f[1]} ){

				print O "$f[1]\t$f[2]\t$f[4]\t$f[6]\t$f[5]\t$kin_true{$f[2]}{$f[1]}\t$pop{$f[2]}$pop{$f[1]}\n";
	    	}
			#else{  print "Not find paris $f[0]<->$f[1] in $input_lcmlkin\n";}
		}
		else{
			print O "Ind1\tInd2\tNSNP\tKin_est\tIBD0_est\trelationship\tKin_true\tIBD0_true\tclass\n";

		}
	}
	close F;
	close O;
}

sub getFormatKINGrobust(){
	$input=$test;
	open F, "<$input";
	open O, ">$output".".plot.fmt";
	while(<F>){
		chomp;
		$line++;
		@f=split(/\s+|\t/,$_);
		$POP=-1;
		if($line>1){
			if( exists $kin_true{$f[2]}{$f[1]} ){
				$IBD0=2*$f[6]/$f[5];
				print O "$f[1]\t$f[2]\t$f[4]\t$f[7]\t$IBD0\t$kin_true{$f[2]}{$f[1]}\t$pop{$f[2]}$pop{$f[1]}\n";
	    	}
			#else{  print "Not find paris $f[0]<->$f[1] in $input_lcmlkin\n";}
		}
		else{
			print O "Ind1\tInd2\tNSNP\tKin_est\tIBD0_est\trelationship\tKin_true\tIBD0_true\tclass\n";

		}
	}
	close F;
	close O;
}


