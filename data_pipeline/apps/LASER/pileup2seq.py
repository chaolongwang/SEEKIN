#!/usr/bin/env python -B
import sys, os, re
try:
    from XiaoweiLib import myopen
except:
    sys.path = [re.sub(r'^/home/zhanxw/', '/net/fantasia/home/zhanxw/', x) for x in sys.path]
    sys.path.append('/net/nfsb/fantasia/home/zhanxw/mylib/Python/')
    sys.path.append(os.path.abspath('.'))
    from XiaoweiLib import myopen

def banner():
    b = """
    ==================================================================
    ====         pileup2seq: Convert pileup to seq format         ====
    ====        Version 1.07, last updated on Feb 22, 2016        ====
    ====      Joint efforts of Xiaowei Zhan and Chaolong Wang     ====
    ====           Feedbacks/comments: zhanxw@gmail.com           ====
    ====         (C) 2013-2016 Xiaowei Zhan, GNU GPL v3.0         ====
    ==================================================================
    """
    print b
    
def usage():
    print("%s [-b bedFileToExcludeRegion] [-i idFile] [-f referenceFile] -m siteFile -o outputPrefix pileupFiles...." % sys.argv[0] )
    print("siteFile: format: chrom, pos, id, ref, alt")
    print("bedfile: format:chrom, beg, end, rsNumber")
    print("idFile: pileupNumber, GWAS_ID - optional")

# return counts (ref, alt)
def count(ref, reads):
    if len(ref) != 1:
        return (0, 0)
    # stripping ^ and $
    reads = re.sub(r'\^.', '', reads)
    reads = reads.replace('$', '')
    reads = reads.replace('$', '')
    
    # stripping insertion and deletion
    while reads.find('+') > 0:
        pos = reads.find('+')
        end = pos + 1
        while reads[end].isdigit():
            end += 1
        reads = reads[:pos] + reads[end:]
    while reads.find('-') > 0:
        pos = reads.find('-')
        end = pos + 1
        while reads[end].isdigit():
            end += 1
        reads = reads[:pos] + reads[end:]
    reads = reads.replace('*','')

    # count ref and alt 
    r, a = 0, 0
    bases = set(list('acgtnACGTN'))
    for i in reads:
        if i == '.' or i == ',':
            r += 1
            continue
        if i in bases:
            a += 1
            continue
        print >> sys.stderr, "Unrecognized base", reads
        break
    return (r, a)

def printHeader(colDict, mapRef, fout):
    row1= ['.' for i in xrange(6)]
    row2= ['.' for i in xrange(6)]
    for pos, rsId in colDict.iteritems():
        ref = mapRef[rsId]
        row1.append(rsId)
        row2.append(ref)
    fout.write('\t'.join(row1))
    fout.write('\n')
    fout.write('\t'.join(row2))
    fout.write('\n')       

# @param qual is a string of character
# character -> error rate
# e.g. phred2qual['!'] = 1
phred2qual = None
from math import log10
def calculateMeanQuality(qual):
    global phred2qual
    if not phred2qual:
        ##logger.info("Update quality table (phred2qual)")
        phred2qual = {}
        for i in xrange(33, 256):
            phred2qual[chr(i)] = 10. ** (- (i - 33.) / 10)
    totalQual = 0.
    totalCount = 0
    for q in qual:
        totalQual +=phred2qual[q]
        totalCount += 1
    if totalCount == 0:
        return 0  ## quality 0
    return round(- 10.*log10(totalQual / totalCount))

def checkReferenceAllels(refFile, siteFile):
    from XiaoweiLib import GenomeSequence
    gs = GenomeSequence()

    chroms = gs.open(refFile)
    logger.info("Load reference file: " + refFile + " (%d chromosomes)" % chroms )

    msg = "Error: Found mismatched reference alleles at %s:%d. [%s in " + siteFile + " vs. %s in " + refFile + "]"
    mismatchedRefAllele = 0
    for ln in myopen(siteFile):
        fd = ln.strip().split()
        if fd[0].lower().startswith('chr'): continue
        chrom, pos, rsid, ref, alt = fd[:5]
        pos = int(pos)
        trueRef = gs.getBase1(chrom, pos)

        if trueRef.lower() != ref.lower():
            mismatchedRefAllele += 1
            logger.warn(msg % (chrom, pos, ref, trueRef))

    if mismatchedRefAllele == 0:
        logger.info("No mismatched reference alleles detected.")
    else:
        logger.info("Detected %d mismatched reference alleles. Please fix and rerun this script." % mismatchedRefAllele )
        sys.exit(1)
    return 0

if __name__ == '__main__':
    banner()
    try:
        import getopt
        optlist, args = getopt.getopt(sys.argv[1:], 'b:i:m:o:f:')
        optlist = dict(optlist)
        arg_bedFile = optlist.get('-b', None)
        arg_idFile = optlist.get('-i', None)
        arg_refFile = optlist.get('-f', None)
        arg_siteFile = optlist['-m']
        arg_outPrefix = optlist['-o']
        fns = args
    except:
        usage()
        raise
        sys.exit(1)

    ## set-up logging
    import logging
    logger = logging.getLogger('pileup2seq')

    fh = logging.FileHandler(arg_outPrefix + '.log')
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s', logger.setLevel(logging.DEBUG))
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    from time import ctime
    logger.info("Started time: " + ctime())
    
    from XiaoweiLib import BedFile
    bedFile = BedFile()
    if arg_bedFile != None:
        logger.info("Load bed lines and unique regions: " +  str(bedFile.open(arg_bedFile, trimChrPrefix = True)))
    else:
        logger.info("Skip loading bed file, all loci will be processed.")

    mapContent = [i.strip().split() for i in open(arg_siteFile).xreadlines()]
    ## site file has header, need to check it
    if len(mapContent) < 1:
        print >>sys.stderr, "Site file is too short or incorrect!"
        sys.exit(1)
    if mapContent[0][:5] != ['CHR', 'POS', 'ID', 'REF', 'ALT'] and mapContent[0][:5] != ['CHROM', 'POS', 'ID', 'REF', 'ALT']:
        print >>sys.stderr, "Site file does not have proper header!!! We try to continue though..."
        print mapContent[0]
    else:
        mapContent = mapContent[1:]
        
    if arg_refFile != None:
        checkReferenceAllels(arg_refFile, arg_siteFile)
    else:
        logger.info("Skip loading reference file, and no reference alleles will be checked")
        
    pos = ['%s:%s' % ( i[0].replace('chr',''), i[1]) for i in mapContent]
    rsId = [ i[2] for i in mapContent]
    # print >> sys.stderr, "Open bed file", bedFile
    # bedContent = [i.split() for i in open(bedFile).xreadlines()]
    #pos = [ '%s:%s' % (i[0].replace('chr',''), i[2]) for i in bedContent]
    #rsId = [ i[3] for i in bedContent]
    # below OrderedDict only in Python 2.7
    from ordereddict import OrderedDict
    colDict = OrderedDict(zip(pos, rsId))  # pos -> rsId e.g. 1:123->rs99
    if len(pos) != len(colDict):
        logger.warn("Total [ %d ] duplicated markers are found in site file [ %s ]" % (len(pos) - len(colDict), arg_siteFile))
        ##print >> sys.stderr, "Duplicated site: ",
        
        ##sys.exit(1)
        dupPos = set()
        processedPos = set()
        for p in pos:
            if p not in processedPos:
                processedPos.add(p)
            else:
                dupPos.add(p)
        if len(dupPos) > 0 :
            logger.warn("Duplicated site: %s " % ",".join([i for i in dupPos]) )
        
    excludePos = [bedFile.contain( i[0].replace('chr',''), i[1]) for i in mapContent]
    excludePos = set([pos[idx] for idx, v in enumerate(excludePos) if v])
    logger.info("Excluding %d on-target markers" % len(excludePos))
    logger.info("%d markers loaded" % len(colDict))
    logger.info("%d pileup files" % len(fns))

    #mapRef = dict([i.strip().split() for i in open(mapRefFile).xreadlines()]) # rsId -> ref
    #mapRef = dict([ ])
    #print >> sys.stderr, "%d refbases loaded" % len(mapRef)

    if arg_idFile != None:
        idContent = []
        pileupId = {} 
        for i in open(arg_idFile).xreadlines():
            i = i.strip().split()
            if i[0].find('.') >= 0:  ## since only before dot part inferred from base file name will be used as key
                logger.error("First column in id file contain '.' as in " + i[0] + ", only before dot parts will be used.")
                i[0] = i[0].split('.')[0]
            idContent.append(i) # pileupId -> GWAS_ID
            val = i[1:]
            if len(val) > 2:
                val = val[:2]
            elif len(val) == 1:
                val = [val[0], val[0]]
            elif len(val) == 2:
                pass
            else:
                logger.error("Wrong line in idFile: "+i)
                continue
            if i[0] in pileupId:
                logger.info("Skipped: Duplicated id in idFile: " + i[0])
                continue
            assert(len(val) == 2), i
            pileupId[i[0]] = val
        ## print pileupId
        logger.info("%d sample id loaded" % len(pileupId))
        
    # covFile = open(outPrefix + '.coverage', 'w')
    # refCountFile = open(outPrefix + '.refCount', 'w')
    # printHeader(colDict, mapRef, covFile)
    # printHeader(colDict, mapRef, refCountFile)
    #logFile = open(arg_outPrefix + '.log', 'w')
    from collections import Counter
    seqFile  = open(arg_outPrefix + '.seq', 'w')
    for fn in fns:
        res = {}
        for ln in myopen(fn):  # loop each pileup file
            fd = ln.strip().split()
            if len(fd) == 4:
                ## after samtools 1.0.18, trancated pileup lines will be outputted
                ## and we will need to ignore them
                chrom, pos, ref, depth = fd
                refCount, altCount, qual = 0, 0, 0
            elif len(fd) == 6:
                chrom, pos, ref, depth, reads, quals = ln.strip().split()

                try:
                    refCount, altCount = count(ref, reads)
                    qual = calculateMeanQuality(quals)
                except:
                    logger.warn("Cannot parse pileup data, entering debug mode ...")
                    import pdb
                    pdb.set_trace()

            else:
                logger.warn("File [ %s ] is empty or not valid, check this line: [ %s ]" % (fn, ln))
                break
            #print 'reads = ', reads,
            #print 'refCnt = ', refCount, 'altCnt=', altCount
            key = '%s:%s' % (chrom.replace('chr', ''), pos )
            if key not in colDict:
                continue
            res[key] = (refCount, altCount, qual)

        # outputs            
        # write out IDs
        if arg_idFile != None:
            key = os.path.basename(fn).split('.')[0]
            gwasId = pileupId.get(key, [key, key]) # value is a list of >= 2 elements
            assert(len(gwasId) == 2), gwasId
            if key in pileupId:
                seqFile.write('\t'.join(gwasId))
            else:
                logger.info("Not update ID: id file does not have correct entry for " + str(key))
                seqFile.write('\t'.join(gwasId))
        else:
            gwasId = os.path.basename(fn).replace('.pileup', '')
            seqFile.write('\t'.join([gwasId, gwasId]))

        # write out ref, alt, qual
        totalDepth = 0
        siteExcluded = 0
        siteHasPileup = 0
        siteNoPileup = 0
        for k, v in colDict.iteritems():
            if k in excludePos:
                seqFile.write('\t0 0 0')
                logger.debug('Exclude\t%s\n' % k)
                siteExcluded += 1
                continue
            if k in res:
                ref, alt, qual = res[k]
                seqFile.write('\t%d %d %d' % (ref+alt, ref, qual))
                siteHasPileup += 1
                totalDepth += (ref+alt)
            else:
                #"no pile up at position"
                seqFile.write('\t0 0 0') # last 0 mean quality zero
                siteNoPileup += 1
        seqFile.write('\n')
        avgDepth = 1.0 * totalDepth / (1e-10 + len(colDict))
        ptgSitePileup = 100.0 * siteHasPileup / (1e-10 + siteHasPileup + siteNoPileup)
        logger.info("%s, avgDepth=%.4f, ptgSiteHasPileup=%.4f%%" % (fn, avgDepth, ptgSitePileup))
    seqFile.close()
    logger.info("Sequence file is generated: " + arg_outPrefix + '.seq')
    
    siteFile  = open(arg_outPrefix + '.site', 'w')
    siteFile.write('\t'.join(['CHR', 'POS', 'ID', 'REF', 'ALT']))
    siteFile.write('\n')
    for fd in mapContent:
        chrom, pos, rsid, ref, alt = fd
        siteFile.write("%s\t%s\t%s\t%s\t%s\n" % (chrom, pos, rsid, ref, alt))
    siteFile.close()
    logger.info("Site file is generated: " + arg_outPrefix + '.site')
    
    ## logFile.close()
    logger.info("Pileup2seq finished at " + ctime())
