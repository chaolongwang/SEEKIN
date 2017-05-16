# SEEKIN software

### SEEKIN:SEequence-based Estimation of KINship.

#### Author: Jinzhuang Dou  <douj@gis.a-star.edu.sg>, Chaolong Wang <wangcl@gis.a-star.edu.sg>

#### License: GNU General Public License v3.0 (GPLv3)
---

## 1 Description
SEEKIN is a software program for estimating pairwise kinship coefficients using either shallow sequencing data or genotyping data. The method was initially developed for shallow sequencing data, such as off-target sequencing data from target sequencing experiments (typically ~0.1-1X). SEEKIN, together with the [LASER](http://csg.sph.umich.edu/chaolong/LASER/) software [1,2] for ancestry estimation, enables control of population structure and cryptic relatedness in target/exome sequencing studies.
 
To address the missing data and genotype uncertainty issues that are intrinsic to shallow sequencing data, we use the LD-based genotype calling algorithm in [BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html) [3] to process sequencing data and develop kinship estimators that explicitly model the genotype uncertainty via the Rsq metric output by [BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html).

SEEKIN includes a homogeneous estimator for samples from the sample population (SEEKIN-hom) and a heterogeneous estimator for samples with population structure and admixture (SEEKIN-het). For SEEKIN-het, we use the [LASER](http://csg.sph.umich.edu/chaolong/LASER/) program to estimate individual ancestry for each study individual, and then derive individual-specific allele frequencies for kinship estimation in a way similar to existing methods [4,5].

Our SEEKIN estimators are also applicable to high-quality genotyping data, such as those from array-genotyping or deep whole-genome sequencing studies. In this special case of no genotype uncertainty, SEEKIN performs similarly to existing kinship estimation programs. 

In terms of implementation, SEEKIN utilizes a "single producer/consumer" design for parallel computation and takes the standard gzipped VCF file as the input.  SEEKIN is both computationally and memory efficient and is scalable to estimate pairwise kinship between 10,000s of individuals on a typical computing cluster.

If you have any bug reports or questions, please send an email to Jinzhuang Dou at <douj@gis.a-star.edu.sg> or Chaolong Wang at <wangcl@gis.a-star.edu.sg>.

## 2 Citation for SEEKIN 

Details of our SEEKIN method can be found at:  
Dou, J. et al. Estimation of kinship coefficients using sparse sequencing data. (Manuscript submitted)


## 3 Dependencies
* gcc >= 4.9
* [OpenBLAS](http://www.openblas.net/)
* [Armadillo](http://arma.sourceforge.net/)

## 4 Download and install

You can download SEEKIN by typing the following command in a terminal.

`git clone https://github.com/jinzhuangdou/SEEKIN.git` 

This command will create a folder named "SEEKIN" in the current directory. The downloaded package contains a statically linked binary executable seekin (in the `bin/` folder), which was pre-compiled and tested in 64-bit Linux machine.  

If you want to compile your own version, please enter the `src/` folder, change the library paths in the beginning of the Makefile accordingly, and type the following command to compile

`cd SEEKIN/src && make`

## 5 Usage 
You can type the following command to get the list of help option.

`seekin -h`  

SEEKIN provides three modules 

* **modelAF** for calculating the PC-related regression coefficients of reference samples;
* **getAF** for estimating the individual allele frequencies of study samples; 
* **kinship** for estimating kinship coefficients for samples from either homogenous or heterogenous samples.  

To get the detailed meaning of option for one module (for example `kinship`), you can type: 

`seekin kinship -h`  

## 6 Examples

Here we provide demo of SEEKIN based on data provided in the `example` folder, which include:

* **Study.10K.vcf.gz**   This file includes genotypes of ~10,000 randomly selected SNP for 50 Chinese and 50 Malays. The genotypes were called from off-target data in a WES study using [BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html) with a reference panel from the 1000 Genome Project. This example dataset is highly noisy and is only for illustration purpose. 
* **SGVP.12K.vcf.gz**   This file includes array genotyping data on ~12,000 SNPs, including ~9,000 SNPs overlapping with the Study.10K.vcf.gz file. This dataset is a subset of the data from the Singapore Genome Variation Project  [SGVP](http://phg.nus.edu.sg/#) [6]. This dataset is intended to be the ancestry reference panel for estimating individual-specific allele frequencies for SEEKIN-het. 
* **SGVP.RefPC.coord**   This file contains PCA coordinates of the top 10 PCs for the reference individuals in SGVP.12K.vcf.gz. This file was generated by the [LASER](http://csg.sph.umich.edu/chaolong/LASER/) software. 
* **Study.onSGVP.PC.coord**   This file contains the top 2 PCs for the study individuals in the SGVP reference ancestry space (`SGVP.RefPC.coord`). This file can be prepared using the [LASER](http://csg.sph.umich.edu/chaolong/LASER/) software with either sparse sequence reads or array genotyping data. 

  
#### 6.1 SEEKIN-hom: kinship estimation for homogenous samples

When assuming no population structure, we only need the genotype file of study samples (`Study.10K.vcf.gz`) to estimate kinship using SEEKIN.   

 ```
  ./seekin kinship -i ./Study.10K.vcf.gz  -r 0.3  -m 0.05   -d DS  -p homo  -l 2000  -t 3  -w  1 -o Study.homo
 ``` 
  
The meaning of all command line options are listed below:

```
Usage: seekin kinship [options]

Options:
      -i  File name of genotypes or dosages of study individuals (gzipped VCF). [no default]
      -f  File name of allele frequencies for study individuals (gzipped VCF). [default "indivAF.vcf.gz"]
      -r  Exclude SNPs with Rsq lower than the -r value. [default  0.5]  
      -m  Exclude SNPs with MAF lower than the -m value. [default 0.05]
      -d  Data used for kinship estimation. GT: genotype values in the GT field; 
          DS: dosage values in the DS field. [default "DS"]
      -p  Kinship estimator. hom: homogeneous estimator; het: heterogeneous estimator. [default "het"]
      -w  Weighting scheme to combine estimates from multiple SNPs. 
          1: MAF(1-MAF)Rsq^2; 2: Rsq^2. [default 1]  
      -l  Number of SNPs in each computing block. [default 10,000]  
      -t  Number of threads for parallel computation. [default 10] 
      -o  Prefix of output file names. 

Note: the option -f is NOT effective when -p is set to "hom". 
Note: when the input file (-i) does not contain Rsq values (e.g. array genotyping data), 
      -d will be set to "GT" and Rsq will be treated as 1.
``` 

The output includes 5 files with prefixes `Study.10K.homo` specified by `-o` flag. 

*  _.log

This is the log file to monitor and record the progress for each module when running SEEKIN.

*  _.kin 

This file provides the kinship estimation for all pairs of individuals.The columns are individual IDs for the first and second individual of pair, number of SNPs used, and the estimated kinship coefficient, respectively. Example: 

  ```
Ind1    Ind2    NSNP    kinship
ME-1M6WZX2      ME-1MYW4LL      8818    0.0051
ME-1M6WZX2      ME-245KCL5      8818    -0.0088
ME-1M6WZX2      ME-2MZQ96C      8818    -0.0006
ME-1M6WZX2      ME-2TGPSG7      8818    -0.0120
ME-1M6WZX2      ME-5P732EC      8818    0.0082
ME-1M6WZX2      ME-667WJSP      8818    -0.0041
ME-1M6WZX2      ME-777J4KU      8818    -0.0023
ME-1M6WZX2      ME-7MLVX2W      8818    -0.0080
ME-1M6WZX2      ME-937Q885      8818    -0.0072
ME-1M6WZX2      ME-BBGMMRQ      8818    0.0070   
  ```
  
*  _.inbreed 

This file provides the estimated inbreeding coefficient for each individual. The first column is individual ID and the second column is the estimated inbreeding coefficient. Example:
  ```
IID     INBREEDCOEF
ME-1M6WZX2      -0.0342
ME-1MYW4LL      -0.0191
ME-245KCL5      -0.0086
ME-2MZQ96C      -0.0150
ME-2TGPSG7      -0.0060
ME-5P732EC      -0.0230
ME-667WJSP      -0.0990
ME-777J4KU      -0.0297
ME-7MLVX2W      -0.0228
ME-937Q885      -0.0807
ME-BBGMMRQ      -0.0219
  ```
  
*  _.matrix and _.matrixID 

The `_.matrix` file contains an N × N matrix of 2$\phi$, where $\phi$ is the estimated kinship matrix of N study individuals. The `.matrixID` file contains the individual IDs corresponding to the order in the `_.matrix` file (also the same order as in the VCF file).


#### 6.2 SEEKIN-het: kinship estimation for heterogenous samples with population structure and admixture

To estimate kinship coefficients for samples with population structure and admixture, we need to first estimate the ancestry background of each individual and derive individual-specific allele frequencies. We use [LASER](http://csg.sph.umich.edu/chaolong/LASER/) to estimate individual ancestry background, which applicable to both shallow sequencing data and genotyping data under a unified framework given an external ancestry reference panel [1,2].  Details of [LASER](http://csg.sph.umich.edu/chaolong/LASER/) can be found on the [LASER](http://csg.sph.umich.edu/chaolong/LASER/). Here, we assume the ancestry coordinates for both the reference individuals and the study individuals have been properly generated by [LASER](http://csg.sph.umich.edu/chaolong/LASER/).

* 6.2.1 Model allele frequencies 

We model allele frequencies as linear functions of PCs based on the ancestry reference panel. The regression coefficients can be estimated using the ```modelAF``` module in SEEKIN:

  ```
  seekin modelAF -i SGVP.12K.vcf.gz  -c SGVP.RefPC.coord -k 2 -o SGVP.beta
  ```
  
This command will generate a file, named SGVP.beta (specified by `–o`), which contains the intercepts and regression coefficients of first two PCs (specified by `-k`) based on genotypes (specified by `-i`) and PC coordinates (specified by `-c`) of the reference individuals. In this file, the first five columns are chromosome, position, reference allele, alternative allele, and allele frequencies of the alternative allele in the reference dataset, respectively. The remaining columns are the estimated intercept and regression coefficients for each PC. This file is tab-delimited. Example:  

  ```
CHROM   POS     REF     ALT     AF      beta0   beta1   beta2
1       740857  G       A       0.013   0.026418      -4.93102e-05    -0.000123628
1       1017197 C       T       0.711   1.421560      -0.000896592     0.000878442
1       1241964 G       A       0.500   0.999935      -0.00123432      0.000292638
1       1462766 C       T       0.351   0.701477      -0.00028479      8.94892e-05
1       1715011 C       A       0.118   0.236869       2.67180000     -0.000551193
1       1979724 C       A       0.226   0.451491      -0.000218917    -0.000157946
1       2251160 T       C       0.397   0.794793       0.000322257    -7.75743e-05
1       2505713 C       T       0.093   0.186592       0.000425398    -0.000163982
1       2749716 G       T       0.308   0.615128      -5.39846e-05    -0.000770491
  ```
  
* 6.2.2 Estimate individual-specific allele frequencies  

The individual specific allele frequencies of study individuals can be generated by the following command:

  ```
  seekin getAF -i Study.onSGVP.PC.coord -b  SGVP.beta  -k 2  -o Study.10K.indvAF.vcf
  ```
This command estimates individual-specific allele frequencies based on the top 2 PCs (specified by `–k`) in the coordinate file (specified by `-i`) and the regression coefficients (specified by `–b`). The outputs are stored in a gzipped VCF file (specified by` –o`). Example:

  ```
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequencies averaged across all individuals">
##FORMAT=<ID=AF1,Number=A,Type=Float,Description="Estimated individual-specific allele frequencies">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ME-UGR9DZ9      ME-RH96Y5Y      ME-WUAGGVA      ME-6Q8D1FD
1       740857  .       G       A       .       .       AF=0.0163       AF1     0.0247  0.0065  0.0269  0.0063  0.0057  0.0160  0.0271
1       1017197 .       C       T       .       .       AF=0.7902       AF1     0.7653  0.8410  0.7530  0.8531  0.8508  0.7571  0.7492
1       1241964 .       G       A       .       .       AF=0.6022       AF1     0.6205  0.6051  0.6188  0.6171  0.6114  0.5651  0.6158
1       1462766 .       C       T       .       .       AF=0.3745       AF1     0.3775  0.3768  0.3767  0.3796  0.3784  0.3657  0.3760
1       1715011 .       C       A       .       .       AF=0.1120       AF1     0.1428  0.0721  0.1519  0.0690  0.0677  0.1179  0.1534
1       1979724 .       C       A       .       .       AF=0.2423       AF1     0.2575  0.2274  0.2607  0.2285  0.2269  0.2376  0.2607
1       2251160 .       T       C       .       .       AF=0.3707       AF1     0.3660  0.3699  0.3664  0.3667  0.3682  0.3804  0.3672
1       2505713 .       C       T       .       .       AF=0.0576       AF1     0.0549  0.0520  0.0565  0.0475  0.0493  0.0709  0.0577
1       2749716 .       G       T       .       .       AF=0.3060       AF1     0.3517  0.2488  0.3646  0.2453  0.2430  0.3116  0.3665
  ```

6.2.3 Estimate kinship coefficients for heterogeneous samples 

With the estimated individual-specific allele frequencies, `Stdudy.10k.indvAF.vcf.gz`, we can run the kinship module to estimate kinship coefficients using the following command:

  ```
  tabix -p vcf ./Study.10K.indvAF.vcf.gz
  
  seekin kinship -i ./Study.10K.vcf.gz  -f  ./Study.10K.indvAF.vcf.gz  -r 0.3  -m 0.05   -d DS  -p het -l 2000  -t 3 -w 1  -o Study.10K.het
  ```
Note that it is important to set `–p` to 'het' so that the SEEKIN-het estimator will be used and the `–f` option will be effective to take the individual-specific allele frequencies as input.
The output files have the same format as those described in section 6.1.



## 6 Reference

1.  Wang, C. et al. Ancestry estimation and control of population stratification for sequence-based association studies. Nat Genet 46, 409-15 (2014)
2.  Wang, C. et al. Improved ancestry estimation for both genotyping and sequencing data using projection procrustes analysis and genotype imputation. Am J Hum Genet 96, 926-937 (2015).
3.  Browning, B.L. & Browning, S.R. A unified approach to genotype imputation and haplotype-phase inference for large data sets of trios and unrelated individuals. Am J Hum Genet 84, 210-23 (2009).
4.  Thornton, T. et al. Estimating kinship in admixed populations. Am J Hum Genet 91, 122-138 (2012).
5.  Conomos, M.P. et al. Model-free estimation of recent genetic relatedness. Am J Hum Genet 98, 127-148 (2016).
6. Teo, Y.Y et al. Singapore Genome Variation Project: a haplotype map of three Southeast Asian populations. Genome Res 19, 2154-2162 (2009).
