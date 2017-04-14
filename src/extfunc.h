#include <unistd.h>
#include <iostream>
#include <cstdlib>
#include <pthread.h>
#include <unistd.h>
#include <cblas.h>
#include <cstdio>
#include <stdio.h>
#include <iomanip> 
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include <sstream>
#include "VcfFileReader.h"
#include "VcfGenotypeField.h"
#include "VcfFileWriter.h"
#include "zlib.h"
#define  ARMA_DONT_USE_WRAPPER    
#include "armadillo"
#include "time.h"
#include "concurrentqueue.h"
#include <pthread.h>
#include <omp.h>
using namespace std;
using namespace arma;
struct readThreadParams {
    fmat G;
    fvec r;
    fmat P;
};



//   


//extern int modelAF ( int arc, char ** argv );
//extern int getAF (int arc, char ** argv);
//extern int kinship(int argc, char ** argv);
