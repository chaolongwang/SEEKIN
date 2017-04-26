
####################################
CC=/opt/gcc-4.9.3/bin/g++ 
## Change $(CC)

ArmadilloPath = /mnt/projects/wangcl/ancestry/usr
## Change $(ArmadilloPath) if Armadillo were installed in a different directory

OpenblasPath =  /mnt/projects/wangcl/ancestry/usr/opt/OpenBLAS
## Change $(OpenblasPath) if OpenblasPath were installed in a different directory
##################################



GCCVERSIONMAJOR := $(shell expr `$(CC) -dumpversion | cut -f1 -d.` \>= 4)
GCCVERSIONMINOR := $(shell expr `$(CC) -dumpversion | cut -f2 -d.` \>= 9)
CFLAGS=         -g -O0 -std=c++11 -pipe 
        
OBJS=		getAF.o  modelAF.o  kinship.o seekin.o 
PROG=       seekin
INCLUDES=   -I ../lib/libStatGen/include  -I $(ArmadilloPath)/include  -I $(OpenblasPath)/include
SUBDIRS=    .
LIBPATH=    -L/lib64 -L/usr/lib64 -L$(ArmadilloPath)/lib64
LIBS=       -pthread ../lib/libStatGen.a  $(OpenblasPath)/lib/libopenblas.a -lgfortran -lz -lm   



BIT_ERR = 0
ifeq (,$(findstring $(shell uname -m), x86_64 ppc64 ia64))
BIT_ERR = 1
endif

.SUFFIXES:.cpp .o

.cpp.o:
		@printf "Compiling $<...                             \r"; \
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< || echo "Error in command: $(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $<"


all:	seekin

modules:
		git submodule update --init --recursive 
envTest:
		@test $(BIT_ERR) != 1 || sh -c 'echo "Fatal: 64bit CPU and Operating System required!";false;'
		@test $(GCCVERSIONMAJOR) == 1 || sh -c 'echo "GCC version lower than 4.4.0";false;'
		@test $(GCCVERSIONMINOR) == 1 || sh -c 'echo "GCC version lower than 4.4.0";false;'

seekin:	envTest $(OBJS)
		@printf "Linking...                                                \r"
		@$(CC)   $(CFLAGS)  -static  -o $(PROG) $(OBJS)  $(LIBPATH) $(LIBS)   $(ENTRAFLAGS) -Xlinker -zmuldefs 
		@printf "$(PROG) compilation done.\n";
		@$(CC) ./*.o  $(LIBPATH) $(LIBS) $(EXTRA_FLAGS) -o ../bin/seekin
		@cp $(PROG) ../bin/

clean:
		@rm -fr gmon.out *.o a.out *.exe *.dSYM $(PROG) *~ *.a *.so.* *.so 
		@printf "$(PROG) cleaning done.\n";

install:
		@printf "$(PROG) installed at ../bin/$(PROG)\n"
