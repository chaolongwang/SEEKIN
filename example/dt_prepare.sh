
LASER_PATH="/mnt/projects/douj/ssmp/software/LASER-2.04"




vcftools --gzvcf /mnt/projects/wangcl/ancestry/resource/LASER-resource/SGVP/SGVP.b37.maf0.01.vcf.gz  --chr 22 --recode --recode-INFO-all --stdout | bgzip -c > SGVP_268.chr22.vcf.gz  

vcftools --gzvcf /mnt/projects/wangcl/merck/kinship/MerckTest/KinshipTestSet/admixture/Merck/result_28031016/0.75_homo_noref/methCompare/All.0.5.vcf.gz --chr 22 --keep  sample.list --recode --recode-INFO-all --stdout | bgzip -c > Study.chr22.vcf.gz 

$LASER_PATH/vcf2geno/vcf2geno --inVcf  SGVP_268.chr22.vcf.gz    --out SGVP_268.chr22
$LASER_PATH/vcf2geno/vcf2geno --inVcf  Study.chr22.vcf.gz    --out Study.chr22
$LASER_PATH/laser  -g  SGVP_268.chr22.geno  -o SGVP_268.chr22  -pca 2 -k 10
$LASER_PATH/trace  -g  SGVP_268.chr22.geno  -s Study.chr22.geno -c SGVP_268.chr22.RefPC.coord  -o Study.chr22  -k 2 -K 4


