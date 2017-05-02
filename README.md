# SEEKIN software

### Sequence-based Estimation of Kinship.

#### Author: Jinzhuang Dou  <douj@gis.a-star.edu.sg>, Chaolong Wang <wangcl@gis.a-star.edu.sg>

#### License: GNU General Public License v3.0 (GPLv3)
---

## 1 Description

SEEKIN is a software program for estimating kinship coefficients for samples which are sequenced at low sequencing coverage (typically lower than 1x). The key features of this program include:   								
* Account for the genotype uncertainties by leveraging the haplotype information from study or external samples.  
* Estimate kinship for heterogenous samples with population structure and admixture.    
* Analyze thousands of individuals in saving memory usage and computational time by utilizing the "single producer/consumer" design. 

If you have any bug reports or questions please send an email to Jinzhuang Dou at `douj@gis.a-star.edu.sg` or Chaolong Wang at `wangcl@gis.a-star.edu.sg`.  

* Citation for SEEKIN:  


## 2 Dependencies
* gcc >= 4.9
* OpenBLAS
* Armadillo

## 3 Download and install

`git clone https://github.com/jinzhuangdou/SEEKIN.git` 

The download package contains a standalone (i.e., statically linked) 64-bit Linux executable seekin (in the `bin/`), which has already been tested on Linux server. If you want to compile your own version of SEEKIN, enter the `src/` folder and type `make` to compile the programs. Before that, you will need to change the library paths in the Makefile accordingly.

`cd SEEKIN/src && make`

## 4 Usage 
You can type the following command to get the list of help option.
`seekin –h`  

SEEKIN provides three modules 

* **modelAF** for calculating the PC-related regression coefficients of reference samples;
* **getAF** for estimating the individual allele frequencies of study samples; 
* **kinship** for estimating kinship coefficients for samples from either homogenous or heterogenous samples.  

To get the detailed meaning of option for one module (for example `kinship`), you can type: `seekin kinship –h`  

## 5 Examples

Here we provide example usages based on the data provided in the folder named `example`. This folder includes two genotype files in the standard VCF compressed format and two PCA coordinate file in the plain text format. 

* **Study.chr22.vcf.gz**  This file includes genotypes of chromosome 22 for 10 studied samples which are called from shallow sequencing reads (~0.75x) and then phased using Beagle (V4.0) software [1] with 1000 Genomes Project Phase 3 (1KG3) as the external reference panel. 
* **SGVP_268.chr22.vcf.gz**  This file includes the genotypes of chromosome 22 for 268 reference samples from Singapore Genome Variation Project (SGVP). 
* **SGVP_268.chr22.RefPC.coord** This file contains PCA coordinates for the top 2 PCs of the reference individuals.
* **Study.chr22.ProPC.coord**  This file contains the top 2 PCs calculated by projecting the study samples on the SGVP panel using LASER [2]. 

  
#### 5.1 Kinship estimation for homogenous samples

Only the genotype file of study samples is required when estimating kinship for homogenous samples.   

  ```./seekin kinship -i ./Study.chr22.vcf.gz  -r 0.3  -m 0.05   -d DS  -p homo  -n 2000  -t 3  -w  1 -o Study.chr22.homo``` 
  
It will generate 5 files with prefixes `Study.chr22.homo` specified by `–o` flag in the current directory. 

*  _.log and terminal outputs 

The terminal outputs are used to monitor and record the progress for each module when running SEEKIN. The log file is identical to the terminal outputs. 

*  _.kin 

This file provides the kinship estimation for all pairs of individuals.The columns are individual ID for the first and second individual of pair, number of SNPs used, the estimated kinship coefficient, respectively. One example is as following: 

  ```
  Ind1    Ind2    NSNP    Kinship      
  S1      S2      8592    0.0231     
  S1      S3      8592    0.0370        
  S1      S4      8592    0.0168      
  ```
  
*  _.inbreed 

This file provides the estimation of inbreeding coefficient estimation for each individual. The columns are individual ID and inbreeding coefficient, respectively. One example is as following:
  ```
  Ind1    Inbreed_coef
  S1       0.0296
  S2       0.0228
  S3      -0.0338
  ```
  
*  _.index and _.matrix 

The `_.matrix` file contains an N × N matrix (The variable N here means the sample size of study samples) of estimated kinship coefficients with the corresponding index of each individual shown in `_.index` file.  For example, the kinship coefficient value given in row 2 and column 3 in the `_.matrix` file would correspond to the individuals in the `_.index` file who have indices of 2 and 3, respectively.

#### 5.2 Kinship estimation of samples with population structure and admixture

We first use the SGVP as reference panel to model the PC-related linear regression coefficients based on the ```modelAF``` module:

  ```seekin modelAF –i SGVP_268.chr22.vcf.gz –c SGVP_268.chr22.RefPC.coord -k 2 –o SGVP_268.chr22.beta```
  
You will generate the file which contains the PC-related regression coefficients for reference samples. From the first column to the fifth column are chromosome ID, genome position, reference allele, alternative reference allele, and allele frequencies of non-ref allele, respectively. And the remaining columns are the estimated coefficients for each PC. This file is tab-delimited. An example is as following: 

  ```
  CHROM   POS     REF     ALT      AF       beta0   beta1   
  10     60969    C       A        0.48     0.96    -0.00    
  10     70969    G       A        0.41     0.96    -0.10    
  ```

Then, the individual allele frequencies of study samples could be generated based on the following command:

  ```
  seekin getAF –i Study.chr22.ProPC.coord -b  SGVP_268.chr22.beta  -k 2 –o Study.chr22.indvAF.vcf
  ```

The generated file is the standard VCF in the compressed format as following: 

  ```
  ##fileformat=VCFv4.2
  ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated Allele Frequencies of all samples">
  ##FORMAT=<ID=AF1,Number=A,Type=Float,Description="Estimated individual specific Allele Frequencies">
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ID1   ID2  ID3
  1       11008   .       C       G       .       .       AF=0.0500       AF1     0.0534  0.0455  0.0536
  1       12001   .       A       G       .       .       AF=0.0200       AF1     0.0231  0.4451  0.1537
  1       13102   .       C       T       .       .       AF=0.4500       AF1     0.2514  0.0123  0.0216
  1       14052   .       C       G       .       .       AF=0.6500       AF1     0.0524  0.0252  0.9531
  ```

Finally, we can run the `kinship` module with the above `Stdudy.chr22.indvAF.vcf.gz` as the input, The command is similar with the homogenous setting but requiring the estimated individual allele frequency file specified by the flag `–a` and using the `admix` mode specified by flag `–p`.

  ```
  seekin kinship -i ./Study.chr22.vcf.gz  -a  ./Study.chr22.indvAF.vcf.gz  -r 0.3  -m 0.05   -d DS  -p admix -n 2000  -t 3 -w 1  -o Study.chr22.admix
  ```
  
 The generated output files have the same format with those from the homogenous setting. 


## 6 Reference

1.  Browning, B.L. & Browning, S.R. A unified approach to genotype imputation and haplotype-phase inference for large data sets of trios and unrelated individuals. Am J Hum Genet 84, 210-23 (2009).
2.  Wang, C. et al. Ancestry estimation and control of population stratification for sequence-based association studies. Nat Genet 46, 409-15 (2014).



