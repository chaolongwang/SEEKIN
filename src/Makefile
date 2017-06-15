
####################################
CXX=/opt/gcc-4.9.3/bin/g++ 
## Change $(CXX)

ArmadilloPath = /mnt/projects/wangcl/ancestry/usr
## Change $(ArmadilloPath) if Armadillo were installed in a different directory

OpenblasPath =  /mnt/projects/wangcl/ancestry/usr/opt/OpenBLAS
## Change $(OpenblasPath) if OpenblasPath were installed in a different directory

OpenblasLibPath =  /mnt/projects/wangcl/ancestry/usr/opt/OpenBLAS/lib
## Change $(OpenblasLibPath) if OpenblasPath were installed in a different directory
##################################

# GCCVERSIONMAJOR := $(shell expr `$(CXX) -dumpversion | cut -f1-2 -d.` \>= 4.9) 
CFLAGS = -g -O3 -std=c++11 -pipe -static
OBJS = getAF.o  modelAF.o  kinship.o seekin.o 
PROG = seekin
INCLUDES = -I ../lib/libStatGen/include  -I $(ArmadilloPath)/include  -I $(OpenblasPath)/include
SUBDIRS = .
LIBPATH = -L/lib64 -L/usr/lib64 -L$(ArmadilloPath)/lib64
LIBS = -pthread ../lib/libStatGen.a  $(OpenblasLibPath)/libopenblas.a -lgfortran -lz -lm   

BIT_ERR = 0
ifeq (,$(findstring $(shell uname -m), x86_64 ppc64 ia64))
BIT_ERR = 1
endif

.SUFFIXES:.cpp .o

.cpp.o:
		@printf "Compiling $<...                             \r"; \
		$(CXX) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< || echo "Error in command: $(CXX) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $<"

all:	seekin

seekin:	$(OBJS)
		@printf "Linking...                                                \r"
		@$(CXX) $(CFLAGS) -o $(PROG) $(OBJS) $(LIBPATH) $(LIBS) $(ENTRAFLAGS) -Xlinker -zmuldefs 
		@printf "$(PROG) compilation done.\n";
		@$(CXX) ./*.o $(LIBPATH) $(LIBS) $(EXTRA_FLAGS) -o ../bin/seekin
		@cp $(PROG) ../bin/

clean:
		@rm -fr gmon.out *.o a.out *.exe *.dSYM $(PROG) *~ *.a *.so.* *.so 
		@printf "$(PROG) cleaning done.\n";

