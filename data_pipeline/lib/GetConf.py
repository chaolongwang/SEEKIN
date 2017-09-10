#!/usr/bin/env python3
"""Creates a config file describing your samples that can be used as
input for phasing pipeline (-i)

"""

#--- standard library imports
#
import subprocess
import sys
import os
import argparse
import logging
import glob
import csv

#--- third-party imports
#
import yaml

#--- project specisfic imports
#
# add lib dir for this pipeline installation to PYTHONPATH
PIPELINE_PATH = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))

LIB_PATH = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
if LIB_PATH not in sys.path:
    sys.path.insert(0, LIB_PATH)

#from readunits import ReadUnit
#from readunits import create_rg_id_from_ru
#from readunits import key_for_readunit

__author__ = "Jinzhuang Dou"
__email__ = "douj@gis.a-star.edu.sg"
__copyright__ = "2016 Genome Institute of Singapore"
__license__ = "The MIT License (MIT)"


# only dump() and following do not automatically create aliases
yaml.Dumper.ignore_aliases = lambda *args: True


# global logger
#logger = logging.getLogger(__name__)
#handler = logging.StreamHandler()
#handler.setFormatter(logging.Formatter(
#    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
#logger.addHandler(handler)


def main():
    """main function
    """

    parser = argparse.ArgumentParser(description=__doc__)

    # generic args
    parser.add_argument('-c', "--conf", required=True,
                        help="the conf files where the parameters are")
    parser.add_argument('-o', "--yaml", required=True,
                        help="Output config (yaml) file")
    parser.add_argument('-f', '--force-overwrite', action='store_true',
                        help="Force overwriting of existing file")
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="Increase verbosity")
    parser.add_argument('-q', '--quiet', action='count', default=0,
                        help="Decrease verbosity")

    args = parser.parse_args()

 #   logger.setLevel(logging.WARN + 10*args.quiet - 10*args.verbose)


 #   if os.path.exists(args.yaml) and not args.force_overwrite:
 #       logger.fatal("Cowardly refusing to overwrite existing file %s", args.yaml)
 #       sys.exit(1)

    # Output the sample configure file
    if not os.path.exists("./jobfiles"):
        os.mkdir(os.path.dirname("./jobfiles/OK")) 

    if not os.path.exists("./snp"):
        os.mkdir(os.path.dirname("./snp/OK")) 

    if not os.path.exists("./laser"):
        os.mkdir(os.path.dirname("./laser/OK")) 

    if not os.path.exists("./seekin"):
        os.mkdir(os.path.dirname("./seekin/OK"))

    if not os.path.exists("./tmp"):
        os.mkdir(os.path.dirname("./tmp/OK"))

    if not os.path.exists("./logs"):
        os.mkdir(os.path.dirname("./logs/OK"))

    if os.path.exists("./log/"):
        os.rmdir("./log")


    samples = dict()
    sites = dict()
    chrs=dict()
    laserBatch = dict()

    # Get the job file 
    perl_script=os.path.join(PIPELINE_PATH, "scripts/GetGP.pl")
    params= 'perl '+ perl_script+ ' -c ' + args.conf
    os.system(params)

    perl_script=os.path.join(PIPELINE_PATH, "scripts/RunLASER.pl")
    params= 'perl '+ perl_script+ ' -c ' + args.conf
    os.system(params)



    # Input the region file for snp calling/phasing 
    with open("./tmp/region.Lst") as fh:
        for line in fh:
            (c, s, e)=line.split('_')
            chrs.setdefault(c, []).append(line.strip('\n'))

    with open("./tmp/laser.batch.lst") as fh:
        for line in fh:
            c=line.strip('\n')
            laserBatch.setdefault(c, []).append("")

    with open("./laser/ID.index") as fh:
        for line in fh:
            (c, s)=line.strip('\n').split(' ')
            samples.setdefault(s, []).append("")

 
    with open(args.yaml, 'w') as fh:
        yaml.dump(dict(samples=samples), fh, default_flow_style=False)
        yaml.dump(dict(sites=sites), fh, default_flow_style=False)
        yaml.dump(dict(chr=chrs), fh, default_flow_style=False)
        yaml.dump(dict(laserBatch=laserBatch), fh, default_flow_style=False)

if __name__ == "__main__":
    main()
