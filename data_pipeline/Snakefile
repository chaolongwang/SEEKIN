# standard library imports
#
import os
import subprocess
import shutil
import glob
import re

# third party imports
#
from snakemake.utils import report

# project specific imports
# #################################################

LIB_PATH = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(workflow.snakefile)),"..", "lib"))

print(LIB_PATH)
if LIB_PATH not in sys.path:
    sys.path.insert(0, LIB_PATH)


# FIX ME 
# non-login bash
shell.executable("/bin/bash")



runChrs=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
laserBatch=config['laserBatch']
print(laserBatch)
CHRS=config['chr']



rule alls:
	input:
		'seekin/seekin.OK'
		

rule seekin:
	input:
		'laser/laser.OK',
		'snp/varCall.OK'
	output:
		'seekin/seekin.OK'
	threads:
		1
	run:
		shell("bash ./jobfiles/runSEEKIN.sh")
		shell("touch {output}")

rule laser:
	input:
        	expand('laser/getPCA{sample}.OK', sample=laserBatch),
	output:
		'laser/laser.OK'
	run:
		shell("bash ./jobfiles/mergeLaserOutput.sh")
		shell("touch {output}")

		
rule laserPerBatch: 
	input:
		'jobfiles/RunLASER{sample}.sh'
	output:
		'laser/getPCA{sample}.OK'
	threads:
		1
	run:
		shell("bash {input}"),
		shell("touch {output}")


rule varCall:
	input:
		expand('snp/getChrGP{chr}.OK', chr=runChrs)
	output:
		'snp/varCall.OK'
	threads:
		1
	run:
		shell("bash ./jobfiles/merge_chr.sh"),
		shell("touch {output}")


rule varCallPerChr: 
 	input: 
 		lambda  wildcards: expand('snp/getGP{job}.OK', job=CHRS[wildcards.chr])
	output:
		'snp/getChrGP{chr}.OK'
	threads:
		1
	run:
		print({input}),
		shell("bash ./jobfiles/merge_chr{wildcards.chr}.sh"),
		shell("touch {output}")
		

rule varCallPerRegion:
	input:
		'jobfiles/getGP_chr{job}.sh'
	output:
		'snp/getGP{job}.OK'
	threads:
		5
	run:
		print ({input}),
		shell("bash {input} 0"),
		shell("touch {output}")
