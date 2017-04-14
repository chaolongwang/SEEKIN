 # Run seekin for homogenous samples
 # ../bin/seekin kinship -i ./Study.chr22.vcf.gz -r 0.3 -m 0.05 -d DS -p homo -n 2000 -t 3 -w 1 -o Study.chr22.homo

 # Run seekin for samples with population structure and admixture
 ../bin/seekin modelAF -i  SGVP_268.chr22.vcf.gz -c SGVP_268.chr22.RefPC.coord -k 2  -o SGVP_268.chr22.beta  

 ../bin/seekin getAF  -i Study.chr22.ProPC.coord -b  SGVP_268.chr22.beta  -k 2  -o Study.chr22.indvAF.vcf

 tabix -p vcf ./Study.chr22.indvAF.vcf.gz


 ../bin/seekin kinship -i ./Study.chr22.vcf.gz  -a  ./Study.chr22.indvAF.vcf.gz  -r 0.3  -m 0.05   -d DS  -p admix -n 2000  -t 3 -w 1  -o Study.chr22.admix
