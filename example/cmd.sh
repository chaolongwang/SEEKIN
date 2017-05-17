

# Run seekin for homogenous samples

../SEEKIN_DEV1.0/bin/seekin kinship -i ./Study.10K.vcf.gz -r 0.3 -m 0.05 -d DS -p hom -l  2000 -t 3 -w 1 -o Study.hom

# Compare the array-based and sequence-based estimates 

Rscript plot.r Study.hom.kin  Study.array.kin  Study.hom



# Run seekin for samples with population structure and admixture

../bin/seekin modelAF -i  SGVP.12K.vcf.gz -c SGVP.RefPC.coord -k 2  -o SGVP.beta  

../bin/seekin getAF  -i Study.onSGVP.PC.coord -b  SGVP.beta  -k 2  -o Study.10K.indvAF.vcf

tabix  ./Study.10K.indvAF.vcf.gz

../bin/seekin kinship -i ./Study.10K.vcf.gz  -f ./Study.10K.indvAF.vcf.gz  -p het  -l 2000 -t 3 -w 1  -o Study.het

# Compare the array-based and sequence-based estimates

Rscript plot.r Study.het.kin  Study.array.kin  Study.het

