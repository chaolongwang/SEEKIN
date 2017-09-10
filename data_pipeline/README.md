SEEKIN Data Preparation and Analysis Pipeline 
=========================================

Summary
----------------------------

This pipeline includes scripts and tools for processing sequencing data to estimate kinship coefficients using SEEKIN. The pipeline is based on [![Snakemake](https://img.shields.io/badge/snakemake-≥3.7.1-brightgreen.svg?style=flat-square)](http://snakemake.bitbucket.org). We have tested the pipeline works on UGE and PBS Pro cluster environments.

The pipeline includes the following steps:

* **Variant calling**

Starting from BAM files, we calculate genotype likelihoods across a given set of SNPs using `samtools mpileup`, followed by LD-based genotype calling using [BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html).   

* **Ancestry estimation**

For each BAM file, we use [LASER](http://csg.sph.umich.edu/chaolong/LASER/) to estimate individual ancestry background with an external ancestry reference panel. 

* **Kinship estimation**

We use [SEEKIN](https://github.com/chaolongwang/SEEKIN) to estimate pairwise kinship coefficients. 


Get started 
----------------------------------------
To run this pipeline, you should prepare a configuration file. One example can be found in `example/example.conf`. In this file, each line specifies one parameter, followed by the parameter value. 

 * BAM_LIST: a text file listing the bam files to analyze, including the file path. Each BAM file corresponds to one study individual. The bam file must be indexed using `samtools index` or equivalent software tools. See `example/sample.bam.list` for an example.
 * VCF_SITE_FILE: a compressed VCF file containing a list of candidate variant sites to be analyzed. One can use the sites.vcf.gz file from [the 1000 Genomes Project](http://www.internationalgenome.org/) with some QC filtering. See `example/example.sites.vcf.gz` for an example.
 * BEAGLE_REF_LIST:  a text file listing the reference panel for BEAGLE analysis (one VCF file per chromosome). See `example/reference.panel.list` for an example.
 
Other parameters can be easily understood according to the comments. More details can be found in the manuals for [SEEKIN](https://github.com/chaolongwang/SEEKIN), [LASER](http://csg.sph.umich.edu/chaolong/LASER/) and [BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html). 

Please remember to modify the path of software installed in your own machine. 

With the configuration file, you can  generate job files by running the following step
```
 $ python $pipelinePath/lib/GetConf.py -c example.conf -o run.yaml
```
After this step, all the jobs required to be run can be found in the folder `./jobfiles`.

To perform variant calling, run the following command 
```
 $ snakemake -s $pipelinePath/Snakefile  --jobs 100  varCall --rerun-incomplete --timestamp --printshellcmds --stats logs/snakemake.stats --configfile run.yaml --latency-wait 60 --cluster-config cluster.yaml --drmaa " -pe OpenMP {threads} -l mem_free={cluster.mem} -l h_rt={cluster.time} -cwd -v PATH -e logs -o logs -w n" --jobname "SEEKIN.slave.{rulename}.{jobid}.sh" >> logs/snakemake.log 2>&1
```
A genotype file called by BEAGLE will be available at `./snp/Beagle.gp.vcf.gz`.

To perform ancestry estimation, run the following command 
```
 $ snakemake -s $pipelinePath/Snakefile  --jobs 10  laser --rerun-incomplete --timestamp --printshellcmds --stats logs/snakemake.stats --configfile run.yaml --latency-wait 60 --cluster-config cluster.yaml --drmaa " -pe OpenMP {threads} -l mem_free={cluster.mem} -l h_rt={cluster.time} -cwd -v PATH -e logs -o logs -w n" --jobname "SEEKIN.slave.{rulename}.{jobid}.sh" >> logs/snakemake.log 2>&1
```
The ancestry coordinates of study samples by LASER will be available at `./laser/laser.seqPC.coord`.

To perform kinship estimation, run the following command 
```
 $ snakemake -s $pipelinePath/Snakefile  --jobs 1  seekin --rerun-incomplete --timestamp --printshellcmds --stats logs/snakemake.stats --configfile run.yaml --latency-wait 60 --cluster-config cluster.yaml --drmaa " -pe OpenMP {threads} -l mem_free={cluster.mem} -l h_rt={cluster.time} -cwd -v PATH -e logs -o logs -w n" --jobname "SEEKIN.slave.{rulename}.{jobid}.sh" >> logs/snakemake.log 2>&1
```
The kinship estimates and other outputs will be available in the folder `./seekin`.

Questions
---------
For further questions, please contact Jinzhuang Dou (douj@gis.a-star.edu.sg) and Chaolong Wang (wangcl@gis.a-star.edu.sg).
