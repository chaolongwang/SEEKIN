system("zless ./snp/chr1.gp.vcf.gz  | grep CHR > ./tmp/vcf.header");
open F, "<./tmp/vcf.header";
while(<F>){
        chomp;
        @f=split(/\s+|\t/,$_);
	$line=0;
	for($i=9;$i<@f;$i++){
	    $line++;
	    $hash{$line}=@f[$i];
	}	
}
close F;
$line=0;
open F, "<./laser/ID.index";
open O, ">./laser/ID.new.index";
while(<F>){
        chomp;
        @f=split(/\s+|\t/,$_);
        $line++;
	print O "$f[0]\t$hash{$line}\n"; 
}
close F;
close O;
system("mv ./laser/ID.new.index  ./laser/ID.index")
