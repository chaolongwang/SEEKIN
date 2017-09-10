"""library functions for read units

following
http://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
"""

# FIXME move read units into own class
# treated as pure dicts from within snakemake though


#--- standard library imports
#
import logging
from collections import namedtuple
from itertools import zip_longest
import hashlib
import os

#--- third-party imports
#
import yaml

#--- project specific imports
#/


__author__ = "Jinzhuang Dou"
__email__ = "douj@gis.a-star.edu.sg"
__copyright__ = "2016 Genome Institute of Singapore"
__license__ = "The MIT License (MIT)"



ReadUnit = namedtuple('ReadUnit',
                      ['run_id', 'flowcell_id', 'library_id', 'lane_id', 'rg_id', 'fq1', 'fq2'])


# global logger
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)


def gen_rg_lib_id(unit):
    """generate read group lib id from readunit"""
    if unit['library_id']:
        return unit['library_id']
    else:
        return "LIB-DUMMY"


def get_sample_for_unit(unitname, config):
    """FIXME:add-doc
    """
    for samplename, readunits in config["samples"].items():
        if unitname in readunits:
            return samplename
    raise ValueError(unitname)


def gen_rg_pu_id(unit):
    """https://www.biostars.org/p/50349/"""
    if unit['run_id'] and unit['flowcell_id'] and unit['lane_id']:
        return "{}_{}.{}".format(unit['run_id'], unit['flowcell_id'], unit['lane_id'])
    else:
        return "PU-" + unit['rg_id']


def fastqs_from_unit(unit):
    """FIXME:add-doc
    """
    if unit['fq2']:
        return unit['fq1'], unit['fq2']
    else:
        return unit['fq1']


def get_Chr_from_cfgfile(cfgfile, raise_off=False):
    """Parse each ReadUnit in cfgfile and return as list
    """

    with open(cfgfile) as fh_cfg:
        yaml_data = yaml.safe_load(fh_cfg)
    unknown_keys = set(yaml_data.keys()) - set(['samples', 'sites'])
    if unknown_keys:
        logger.critical("Found unexpected keys in %s (only 'samples'"
                        " and 'readunits' allowed): %s", cfgfile, unknown_keys)
        if not raise_off:
            raise ValueError(cfgfile)
    samples, readunits_plain = yaml_data['samples'], yaml_data['readunits']

    #logger.debug("samples: {}".format(samples))
    #logger.debug("readunits_plain keys: {}".format(readunits_plain.keys()))
    for sample_key, sample_rus in samples.items():
        for ru_key in sample_rus:
            if ru_key not in readunits_plain.keys():
                logger.critical("readunit %s of sample %s not found"
                                " in config file", ru_key, sample_key)
                if not raise_off:
                    raise ValueError(cfgfile)

    readunits = dict()# actual namedtuples instead of dict
    for ru_key, ru_plain in readunits_plain.items():
        for f in ['run_id', 'flowcell_id', 'library_id', 'lane_id']:
            if f not in ru_plain:
                logger.fatal("Missing field %s in config file %s", f, cfgfile)
                if not raise_off:
                    raise ValueError(cfgfile)
        run_id = ru_plain.get('run_id')
        flowcell_id = ru_plain.get('flowcell_id')
        library_id = ru_plain.get('library_id')
        lane_id = ru_plain.get('lane_id')
        rg_id = ru_plain.get('rg_id')# allowed to be none or missing
        fq1 = ru_plain.get('fq1')
        fq2 = ru_plain.get('fq2')

        # if we have relative paths, make them abs relative to cfgfile
        if not os.path.isabs(fq1):
            fq1 = os.path.abspath(os.path.join(os.path.dirname(cfgfile), fq1))
        if fq2 and not os.path.isabs(fq2):
            fq2 = os.path.abspath(os.path.join(os.path.dirname(cfgfile), fq2))

        for f in [fq1, fq2]:
            if f and not os.path.exists(f):
                logger.fatal("Non-existing input file %s in config file %s", f, cfgfile)
                if not raise_off:
                    raise ValueError(cfgfile)

        ru = ReadUnit(run_id, flowcell_id, library_id, lane_id, rg_id,
                      fq1, fq2)
        if not rg_id:
            ru = ru._replace(rg_id=create_rg_id_from_ru(ru))
        readunits[ru_key] = dict(ru._asdict())

    return samples, readunits


