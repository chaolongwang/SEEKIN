import os, sys, re
try:
    from subprocess import Popen, PIPE, call
except:
    sys.path = [re.sub(r'^/home/zhanxw/', '/net/dumbo/home/zhanxw/', x) for x in sys.path]
    sys.path.append('/net/dumbo/home/zhanxw/python27/lib/python2.7/')
    from subprocess import Popen, PIPE, call
import gzip
import signal

def myopen(fn):
    import gzip
    f = gzip.open(fn)
    try:
        f.read(2)
        f.close()
        return gzip.open(fn)
    except:
        f.close()
        return open(fn)
    
# return merged results as a new dict
def mergeRegion(region):
    for chrom in region.iterkeys():
        ex = region[chrom]
        ex.sort(key = lambda x: x[0])

        ret = []
        for i, r in enumerate(ex):
            if len(ret) == 0:
                ret.append(r)
                continue
            last_end = ret[-1][1]
            beg = r[0]
            end = r[1]
            if beg < last_end:
                if end > last_end:
                    ret[-1][1] = end
            else:
                ret.append(r)
        region[chrom] = ret
    return region

class BedFile:
    def __init__(self):
        self.data = dict()
    # return len of chroms and total bases
    def open(self, fn, trimChrPrefix = False):
        for ln in myopen(fn):
            fd = ln.strip().split()
            chrom, beg, end = fd[:3]
            if chrom[:3].lower() == 'chr':
                if trimChrPrefix:
                    chrom = chrom[3:]
                else:
                    print >> sys.stderr, "Be cautious of 'chr' as chromosome prefix"
            beg, end = int(beg), int(end)
            if chrom in self.data:
                self.data[chrom].append( [beg, end] )
            else:
                self.data[chrom] = [ [beg, end] ]
        self.data = mergeRegion(self.data)
        return len(self.data), sum( (len(self.data[chrom]) for chrom in self.data.iterkeys())) 
    # boundary is left inclusive, right exclusive.
    def contain(self, chrom, pos):
        if chrom not in self.data:
            return False
        try:
            pos = int(pos)
        except:
            return False
        region = self.data[chrom]
        lo = 0
        hi = len(region)
        while lo < hi:
            mid = (lo + hi) /2
            midval = region[mid]
            if pos < midval[0]:
                hi = mid
            elif midval[0] <= pos < midval[1]:
                return True
            elif pos >= midval[1]:
                lo = mid + 1
        return False
    # return a tuple ( a, b, c) where
    # a: whether chrom:pos is contained
    # (b, c): if a == True, get region
    #         if a == False, get distance to left and distance to right
    # e.g.  region: [100, 200), [300,400)
    # return (True, 100, 200) if call getDistance(chrom, 150)
    # return (False, 1, 99) if call getDistance(chrom, 201)
    def getDistance(self, chrom, pos):
        if chrom not in self.data:
            return False
        try:
            pos = int(pos)
        except:
            return (False, None, None)
        region = self.data[chrom]
        lo = 0
        hi = len(region)
        mid = None
        if hi == 0: 
            return (False, None, None)
        while lo < hi:
            mid = (lo + hi) /2
            midval = region[mid]
            if pos < midval[0]:
                hi = mid
            elif midval[0] <= pos < midval[1]:
                return (True, midval[0], midval[1])
            elif pos >= midval[1]:
                lo = mid + 1
        if mid == 0:
            midval = regin[mid]
            
        elif mid == len(region) - 1:
            midval = regin[mid]
        else:
            midval = regin[mid]

        return False

def run(cmd):
    print("= %s" % cmd)
    if (cmd.find(">") <0 ):     # no need to redirect output
        call(cmd.split())
        return
    else:
        p=Popen(cmd, shell=True)
        os.waitpid(p.pid, 0)
        return

def runPool(cmdList, poolSize=4):
    pool = Pool(processes = poolSize)
    pool.map(run, cmdList)

def fastCompressGzip(fn):
    import subprocess
    p = subprocess.Popen(["gzip", "-c", fn], stdout = subprocess.PIPE, close_fds=True, universal_newlines = False)
    return p.stdout

def fastDecompressGzip(fn):
    import subprocess
    p = subprocess.Popen(["zcat", fn], stdout = subprocess.PIPE, close_fds = False, universal_newlines = False)
    return p.stdout
    
"""iterator help you get a list with four elements of any FASTQ record"""
class FastqReader:
    def __init__(self, fileName):
        if (len(fileName)>3 and fileName[-3:]==".gz"):
            self.f = gzip.open(fileName, "rb")
        else:
            self.f = open(fileName)
        signal.signal(signal.SIGPIPE, signal.SIG_DFL) 

    def __iter__(self):
        return self

    def next(self):
        record = [ self.f.readline().strip() for i in range(4)]
        if (record[0] == ''):
            self.f.close()
            raise StopIteration
        return record

class SAMReader:
    def __init__(self, fileName, isBam = False):
        if (isBam == False and fileName[-4:] == ".sam"):
            self.f = open(fileName)
        elif (isBam == True and fileName[-4:] == ".bam"):
            args = "samtools view " + fileName
            print args.split()
            self.f = Popen( args.split(), stdout=PIPE).stdout
        else:
            print("your suffix and does not match value self.isBam")
            sys.exit(1)
        self.header = list() # store SAM header
        self.lineNo = 0 # store current read line number
        self.line = "" # store current read line
        signal.signal(signal.SIGPIPE, signal.SIG_DFL) 
        
    def __iter__(self):
        return self
    def next(self):
        self.line = self.f.readline()
        self.lineNo += 1
        while (self.line != '' and self.line[0] == '@'):
            self.header.append(self.line)
            self.line = self.f.readline()
            self.lineNo += 1

        if (self.line == ''):
            self.f.close()
            raise StopIteration

        fields= self.line.split('\t')
        record = dict()
        record["QNAME"] = fields[0]
        record["FLAG" ] = int(fields[1])
        record["RNAME"] = fields[2]
        record["POS"  ] = int(fields[3])
        record["MAPQ" ] = int(fields[4])
        record["CIGAR"] = fields[5]
        record["MRNM" ] = fields[6]
        record["MPOS" ] = int(fields[7])
        record["ISIZE"] = int(fields[8])
        record["SEQ"  ] = fields[9]
        record["QUAL" ] = fields[10]
        record["TAGS" ] = fields[11:]

# we don't care the optional tags unless necessary
#         if (len(fields) > 11):
#             for i in fields[11:]:
#                 (tag, vtype, value) = i.split(":")
#                 if (vtype=="i"):
#                     record[tag] = int(value)
#                 elif(vtype=="f"):
#                     record[tag] = float(value)
#                 else:
#                     record[tag] = value
        return record

    def dump(self):
        print self.line

# from Dive into Python
def info(object, spacing=10, collapse=1):   
    """Print methods and doc strings.
    
    Takes module, class, list, dictionary, or string."""
    methodList = [method for method in dir(object) if callable(getattr(object, method))]
    processFunc = collapse and (lambda s: " ".join(s.split())) or (lambda s: s)
    print "\n".join(["%s %s" %
                      (method.ljust(spacing),
                       processFunc(str(getattr(object, method).__doc__)))
                     for method in methodList])


#import __main__ # will access the global variables set in __main__
# e.g. Note, you need to put code under __main__ 
#     boolParam = False
#     intParam = 1
#     floatParam = 2.3
#     strParam = "empty"

#     getOptClass = GetOptClass()
#     argumentList = (
#         ('boolParam', ('-b','--bool')),
#         ('intParam',  ('-i', '--integer')),
#         ('floatParam',  ('-f', '--float')),
#         ('strParam', ('-s', '--str'))
#         )
#     getOptClass.parse(argumentList, verbose= False )
#     print 'boolParam = ', boolParam
#     print 'intParam = ', intParam
#     print 'floatParam = ', floatParam
#     print 'strParam = ', strParam
#     print 'rest arguments = ', ",".join(getOptClass.rest)
class GenomeSequence:
    gs = None
    def __init__(self):
        pass
    def isGzFile(self, fn):
        try:
            with open(fn) as f:
                GZMagicNumber = '\x1f\x8b\x08'
                if f.read(3) == GZMagicNumber:
                    return True
        except:
            pass
        return False
    #' Check if it is plain .fa file
    def isPlainFile(self, fn):
        try:
            with open(fn) as f:
                if f.read(1) == ">":
                    return True
        except:
            pass
        return False
    # return number of chromosomes loaded
    def open(self, fn):
        if not os.path.exists(fn):
            print >> sys.stderr, "Reference genome file does not exist: " + fn
            return -1
        
        # check if .fai file exists
        if os.path.exists(fn + '.fai'):
            # check if fn is a gzip file
            if self.isGzFile(fn):
                self.gs = BgzipGenomeSequence()
            elif self.isPlainFile(fn):
                self.gs = PlainIndexedGenomeSeqeunce()
            else:
                print >> sys.stderr, "Unsupported genome sequence type!"
                return -1
        else:
            self.gs = InMemoryGenomeSequence()
        return self.gs.open(fn)
    def read(self, fn):
        return self.gs.open(fn)
    # pos: 0-based index
    def getBase0(self, chrom, pos):
        return self.gs.getBase0(chrom, pos)
    # pos: 1-based index
    def getBase1(self, chrom, pos):
        return self.gs.getBase1(chrom, pos)
        
class InMemoryGenomeSequence:
    gs = dict()
    def __init__(self):
        pass
    def open(self, fn):
        content = [ln.strip() for ln in myopen(fn).readlines() if len(ln.strip()) > 0 ]
        chromIdx = [i for i, ln in enumerate(content) if ln[0] == '>' ]
        chroms = [content[i][1:].split()[0].replace('chr', '') for i in chromIdx]
        seqIdx = [ (chromIdx[i], chromIdx[i+1]) for i in xrange(len(chromIdx) - 1) ]
        seqIdx.append( (chromIdx[-1], len(content) ) )
        seq = [ ''.join(content[ ( i[0] + 1) : i[1]]) for i in seqIdx]
        self.gs = dict(zip(chroms, seq))
        for k, v in self.gs.iteritems():
            print >> sys.stderr, "Chromosome %s loaded with %d bases" % (k, len(v))
        return len(self.gs)
    def read(self, fn):
        return self.open(fn)
    # pos: 0-based index
    def getBase0(self, chrom, pos):
        chrom = chrom.replace('chr','')
        pos = int(pos)
        if chrom not in self.gs:
            return 'N'
        return self.gs[chrom][pos]
    # pos: 1-based index
    def getBase1(self, chrom, pos):
        chrom = chrom.replace('chr','')
        pos = int(pos)
        if chrom not in self.gs:
            return 'N'
        return self.gs[chrom][pos - 1]

class PlainIndexedGenomeSeqeunce:
    handle = None
    index = {}
    def __init__(self):
        pass
    def open(self, fn):
        self.handle = open(fn, 'r')
        # read fai
        for ln in myopen(fn + '.fai'):
            fd = ln.strip().split()
            self.index[fd[0].replace('chr', '')] = [int(i) for i in fd[1:]]
        return len(self.index)            
    def read(self,fn):
        return self.open(fn)
    def getBase0(self, chrom, pos):
        chrom = chrom.replace('chr','')
        if chrom not in self.index:
            return 'N'
        chromLen, fileOffset, nchar, nchar2 = self.index[chrom]
        if pos >= chromLen or pos < 0:
            return 'N'
        a, b = divmod(pos, nchar)
        offset = fileOffset + a * nchar2 + b
        self.handle.seek(offset)
        return self.handle.read(1)

    # pos: 1-based index
    def getBase1(self, chrom, pos):
        return self.getBase0(chrom, pos - 1)

class BgzipGenomeSequence:
    handle = None
    index = {}
    def __init__(self):
        pass
    def open(self, fn):
        try:
            from Bio import bgzf
        except:
            print >> sys.stderr, "Cannot import Bio.bgzf, need to check the installation"
        self.handle = bgzf.BgzfReader(fn)
        # read fai
        for ln in myopen(fn + '.fai'):
            fd = ln.strip().split()
            self.index[fd[0].replace('chr', '')] = [int(i) for i in fd[1:]]
        return len(self.index)
    def read(self,fn):
        return self.open(fn)
    def getBase0(self, chrom, pos):
        chrom = chrom.replace('chr','')
        if chrom not in self.index:
            return 'N'
        chromLen, fileOffset, nchar, nchar2 = self.index[chrom]
        if pos >= chromLen or pos < 0:
            return 'N'
        a, b = divmod(pos, nchar)
        offset = fileOffset + a * nchar2 + b
        self.handle.seek(offset)
        return self.handle.read(1)
    # pos: 1-based index
    def getBase1(self, chrom, pos):
        return self.getBase0(chrom, pos - 1)
    
class GetOptClass:
    rest = []
    def __init__(self):
        pass
    def parse(self, optList, verbose=False):
        # assigne default values should be finished outside of this class
        # preprocess argument
        optDict = dict()
        for i in optList:
            for j in i[1]:
                if len(j) == 0 or j[0] != '-':
                    print "Illegal options: %s" % j
                    sys.exit(1)
                optDict[j] = i[0]

        # parse sys.argv
        index = 1
        while index < len(sys.argv):
            opt = sys.argv[index]
            if opt in optDict:
                # assign the parameter
                varName = optDict[opt]
                # detect if the parameter is bool
                isBoolType = False
                exec('isBoolType = False.__class__ == __main__.%s.__class__' % varName )
                #print isBoolType
                if (isBoolType):
                    try:
                        exec("__main__.%s = True" % varName ) 
                    except:
                        print "Cannot set the boolean variable: %s" % varName
                else:
                    try:
                        arg = sys.argv[index+1]
                    except:
                        print "Not provided argument for option: %s" % varName
                    try:
                        exec("__main__.%s = (__main__.%s.__class__)(%s)" % (varName, varName, repr(arg)))
                    except:
                        print "Cannot set the boolean variable: %s" % varName
                    index += 1
            else:
                self.rest.append(opt)

            index += 1
        
        #check result
        if (verbose == True):
            self.dump(optList, optDict)
        
    def dump(self, optList, optDict):
        print "%s by Xiaowei Zhan" % sys.argv[0]
        print
        print "User Specified Options"
        for i in optList:
            try:
                print (",".join(i[1])).rjust(20),':',
                exec ("print %s " % optDict[i[1][0]])
            except:
                print "Failed to dump ", i
                raise
        print 'rest arguments'.rjust(20), ':', ",".join(self.rest)

# from Dabeaz's great slides
import os 
import fnmatch
def gen_find(filepat,top): 
    for path, dirlist, filelist in os.walk(top):
        for name in fnmatch.filter(filelist,filepat): 
            yield os.path.join(path,name)

import gzip, bz2 
def gen_open(filenames):
    for name in filenames: 
        if name.endswith(".gz"):
            yield gzip.open(name) 
        elif name.endswith(".bz2"):
            yield bz2.BZ2File(name) 
        else:
            yield open(name)

def gen_cat(sources): 
    for s in sources:
        for item in s: 
            yield item

def gen_grep(pat, lines): 
    patc = re.compile(pat) 
    for line in lines:
        if patc.search(line): 
            yield line

