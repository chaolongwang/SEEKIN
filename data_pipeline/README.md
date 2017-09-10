Sequence-based Kinship Estimation Pipeline 
=========================================

Summary
----------------------------

This pipeline includes the software tools for estimating pairwise kinship coefficients starting from the sequencing datasets. It is based on [![Snakemake](https://img.shields.io/badge/snakemake-≥3.7.1-brightgreen.svg?style=flat-square)](http://snakemake.bitbucket.org). Now, the pipeline works on GIS aquila (UGE) or the National Super Computing Center (NSCC) (PBS Pro) environment.
The following steps are performed:

* **Variant calling**

To address the missing data and genotype uncertainty, we call the variants by `samtools mpileup` tools and then process them based on the LD-based genotype calling algorithm in [BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html).   

* **Ancestry estimation**

For each sequenced genome (in BAM), we use [LASER](http://csg.sph.umich.edu/chaolong/LASER/) to estimate individual ancestry background given an external ancestry reference panel. 

* **Kinship estimation**

We use [SEEKIN](https://github.com/chaolongwang/SEEKIN) to estimate pairwise kinship coefficients for homogenous/heterogenous samples. 


Steps to perform kinship estimation 
----------------------------------------
To estimate kinship coefficients using this pipeline, you should prepare for the configuration file. One example can be seen in `example/example.conf`. In this file, each line specifies one parameter, followed by the parameter value. 

 * BAM_LIST: aligned sequenced reads in BAM format. Each BAM file should contain one sample per subject. It also must be indexed using `samtools index` or equivalent software tools. See `example/sample.bam.list` for example.
 * VCF_SITE_FILE: candidate variant sites file in the VCF format. This file includes region in which `samtools mpileup` is generated. This file can include the markers from [the 1000 Genomes Project](http://www.internationalgenome.org/). See `example/example.sites.vcf.gz` for example.
 * BEAGLE_REF_LIST:  external reference panel file list for beagle imputation (one VCF file per chromosome). See `example/reference.panel.list` for example.
 
Other parameters are easily understood according to the comments. More details can be seen in the [SEEKIN](https://github.com/chaolongwang/SEEKIN), [LASER](http://csg.sph.umich.edu/chaolong/LASER/) and [BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html) manuals. Please remember to modify the path of software to specify it installed in your own machine. 

Then, you can  generate the job files by running the following step
```
 $ python  $pipelinePath/lib/GetConf.py  -c example.conf   -o run.yaml
```
After this step, all the jobs required to be run can be seen in the folder `./jobfiles`.

To perform variant calling, run the following step 
```
 $ snakemake -s $pipelinePath/Snakefile  --jobs 100  varCall --rerun-incomplete --timestamp --printshellcmds --stats logs/snakemake.stats --configfile run.yaml --latency-wait 60 --cluster-config cluster.GIS.yaml --drmaa " -pe OpenMP {threads} -l mem_free={cluster.mem} -l h_rt={cluster.time} -cwd -v PATH -e logs -o logs -w n" --jobname "SEEKIN.slave.{rulename}.{jobid}.sh" >> logs/snakemake.log 2>&1
```
The generated genotype file will be available at `./snp/Beagle.gp.vcf.gz`.

To perform ancestry estimation, run the following step 
```
 $ snakemake -s $pipelinePath/Snakefile  --jobs 10  laser --rerun-incomplete --timestamp --printshellcmds --stats logs/snakemake.stats --configfile run.yaml --latency-wait 60 --cluster-config cluster.GIS.yaml --drmaa " -pe OpenMP {threads} -l mem_free={cluster.mem} -l h_rt={cluster.time} -cwd -v PATH -e logs -o logs -w n" --jobname "SEEKIN.slave.{rulename}.{jobid}.sh" >> logs/snakemake.log 2>&1
```
The generated PCA coordinate file of study samples will be available at `./laser/laser.seqPC.coord`.

To perform kinship estimation, run the following step 
```
 $ snakemake -s $pipelinePath/Snakefile  --jobs 1  seekin --rerun-incomplete --timestamp --printshellcmds --stats logs/snakemake.stats --configfile run.yaml --latency-wait 60 --cluster-config cluster.GIS.yaml --drmaa " -pe OpenMP {threads} -l mem_free={cluster.mem} -l h_rt={cluster.time} -cwd -v PATH -e logs -o logs -w n" --jobname "SEEKIN.slave.{rulename}.{jobid}.sh" >> logs/snakemake.log 2>&1
```
The generated output will be available at `./seekin`.

Questions
---------
For further questions, please contact Jinzhuang Dou (douj@gis.a-star.edu.sg) and Chaolong Wang (wangcl@gis.a-star.edu.sg).









 
