/*
 *  Copyright (C) 2015-2015  Genome institute of Singapore (GIS),
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "extfunc.h"

const  string default_str = "NA";
const  int default_int = -999999998;
const  double default_double = -9.99999999;
const  double DELTA=0.001;
static time_t rawtime;
static float runningtime, t1, t2;
static struct tm * timeinfo;
ofstream foutLog;

string INPUT_FILE =default_str;
string OUTPUT_FILE = default_str;
string AF_REF_FILE = default_str;
string AF_INDIV_FILE = "NA";
string POP_STRUCT = default_str;
int WEIGHT=1;
string LOG_FILE="seekin";
string GENE_MODE = "DS"; 
double MAF_FILTER = 0.05;
double R2_FILTER = 0.3;
int BLOCK_SIZE=10000;
int THREAD_NUM=1;

map<int,int> INDEX;
map<int,string> overlapSNPlist; 
map<int,double> overlapRsq; 
map<int,double> overlapMAF; 
map<string,int> geno2dosageMap;
map<string,double> genoCode;


fmat V1;
int overlapSNPcnt=0;
int gp_sampleNum;  
int af_indiv_sampleNum; 
int snp_count=0;
int refAFlag=0; 
int indivAFlag=0; 
int dosageFlag=0;
int rsqFlag=0;
int mode=0;

static bool display_usage();
static void initenv ( int argc, char ** argv );
static vector<string> split(const string& src, string sp);
static string showtime();
static void paraCheck();
static int getVCFsampleNum(string inputfile);
static void checkHeader(int af_indiv_sampleNum, int gp_sampleNum, map<int,int> & INDEX, string gp_file, string af_indiv_file);
static map<string,int> getSNPlist(string inputfile);
static fmat admixEst(string gp_file, string af_indiv_file, int gp_sampleNum, int BLOCK_SIZE, map<string, int> LIST, map <int, int> INDEX);
static fmat homoEst(string gp_file, int gp_sampleNum, int BLOCK_SIZE);
void blockWorkerAdmix(int & BLOCK_SIZE, int & gp_sampleNum, int &numWorking, atomic<bool>& queueFilled,
                    moodycamel::ConcurrentQueue<readThreadParams>& workQueue,
                    moodycamel::ConcurrentQueue<fmat> & resultQueue);

void blockWorkerHomo(int & BLOCK_SIZE, int & gp_sampleNum, int &numWorking,
                    atomic<bool>& queueFilled,
                    moodycamel::ConcurrentQueue<readThreadParams>& workQueue,
                    moodycamel::ConcurrentQueue<fmat> & resultQueue);

static void output(fmat &KIN, string OUTPUT_FILE, int gp_sampleNum, string inputfile, int snp_count);


int kinship(int argc, char ** argv){

    fmat KIN;
    // check the input parameters and set default values. 
    initenv (argc, argv);
    paraCheck();
    openblas_set_num_threads(THREAD_NUM);
    if(OUTPUT_FILE.compare(default_str)==0) OUTPUT_FILE = "seekin";
    if(MAF_FILTER==default_double) MAF_FILTER = 0.05;
    if(R2_FILTER==default_double) R2_FILTER = 0.3;  

    cout<<"["<< showtime() << "] MAF filter " << MAF_FILTER << "\n";
    cout<<"["<< showtime() << "] R2 filter " <<  R2_FILTER << "\n";
    cout<<"["<< showtime() << "] GENOTYPE MODE  " <<  GENE_MODE << "\n";

    foutLog<<"["<< showtime() << "] MAF filter " << MAF_FILTER << "\n";
    foutLog<<"["<< showtime() << "] R2 filter " <<  R2_FILTER << "\n";

    gp_sampleNum=getVCFsampleNum(INPUT_FILE);
    if(POP_STRUCT.compare("homo")==0){
        KIN=homoEst(INPUT_FILE, gp_sampleNum, BLOCK_SIZE);
    }
    else{
        af_indiv_sampleNum=getVCFsampleNum(AF_INDIV_FILE);
        checkHeader(af_indiv_sampleNum, gp_sampleNum, INDEX, INPUT_FILE, AF_INDIV_FILE);
        map<string,int> LIST = getSNPlist (AF_INDIV_FILE);
        KIN=admixEst(INPUT_FILE, AF_INDIV_FILE, gp_sampleNum, BLOCK_SIZE, LIST, INDEX);
    }

    output(KIN, OUTPUT_FILE, gp_sampleNum, INPUT_FILE, overlapSNPcnt);
    //t2 = omp_get_wtime();
    //runningtime =(t2- t1);
    //runningtime=0;
    cout<<"["<< showtime() << "] Finished!\n";
    //cout << "Total CPU time: " << runningtime << " seconds.\n";
    foutLog<<"["<< showtime() << "] Finished!\n";
    //foutLog << "Total CPU time: " << runningtime << " seconds.\n";
    foutLog.close();
}


static fmat homoEst(string gp_file, int gp_sampleNum, int BLOCK_SIZE){

    moodycamel::ConcurrentQueue<fmat> resultQueue;
    moodycamel::ConcurrentQueue <readThreadParams> workQueue;
    vector<thread> workers;
    atomic<bool> queueFilled{false};
    readThreadParams readParams;
    fmat kin=zeros <fmat> (gp_sampleNum, gp_sampleNum);
    fvec r = zeros <fvec> (BLOCK_SIZE);  
    fmat G = zeros <fmat> (BLOCK_SIZE, gp_sampleNum);   

    int numWorking=1;    // Here only one process queue is used
    workers.emplace_back(blockWorkerHomo,ref(BLOCK_SIZE), ref(gp_sampleNum),ref(numWorking),
            ref(queueFilled),ref(workQueue),ref(resultQueue));

    VcfFileReader reader; 
    VcfHeader header;
    VcfRecord  record;
    reader.open(INPUT_FILE.c_str(),header);  
    
    int lineCNT=0;
    int blockReadIsOver=0;
    int matPtr=0;
    bool PASS="TRUE";
    double af_tmp=0;
    double r2_tmp=0;
    cout<<"["<< showtime() << "] Reading the markers in " << INPUT_FILE.c_str() << "\n";
    foutLog<<"["<< showtime() << "] Reading the markers in " << INPUT_FILE.c_str() << "\n";

    while(reader.readRecord(record)){                  
        VcfRecordInfo& info = record.getInfo();
        string chr = record.getChromStr();
        VcfRecordGenotype& geno = record.getGenotypeInfo();
        stringstream pos;
        pos << record.get1BasedPosition();
        string ID = chr+"_"+pos.str();   
        lineCNT++;
        fvec _G = zeros <fvec> (gp_sampleNum);  
        if(GENE_MODE.compare("DS")==0){
            r2_tmp = atof((*(info.getString("DR2"))).c_str());
            af_tmp=0;
            for (int k=0; k<gp_sampleNum; k++){                  
                _G[k]= atof((*geno.getString("DS",k)).c_str());
                af_tmp=af_tmp+_G[k]/2;
            }
            af_tmp=af_tmp/gp_sampleNum;     
        }
        else{ 
            af_tmp=0;
            for (int k=0; k<gp_sampleNum; k++){   
                string tmp= (*geno.getString("GT",k)).c_str();
                _G[k]=genoCode[tmp];
                af_tmp=af_tmp+genoCode[tmp]/2;
            }
            r2_tmp=1;
            af_tmp=af_tmp/gp_sampleNum;
        }
   
        if(r2_tmp >= R2_FILTER & af_tmp >= MAF_FILTER & af_tmp <= (1-MAF_FILTER)) {  //***
            overlapSNPcnt++;
            matPtr=overlapSNPcnt%BLOCK_SIZE;
            if(matPtr==0){
                matPtr=BLOCK_SIZE;
                blockReadIsOver=1;
            }
            overlapSNPlist[matPtr-1]=ID;
            G.row(matPtr-1)=_G.t();
            r[matPtr-1]=r2_tmp;
        }
        if(blockReadIsOver){
            cout<<"["<< showtime() << "] ...Scanning " << lineCNT << " markers and overlapping markers " << overlapSNPcnt<< ".\n";
            foutLog<<"["<< showtime() << "] ...Scanning " << lineCNT << " markers and overlapping markers " << overlapSNPcnt<< ".\n";
            readParams.G=G;
            readParams.r=r;
            workQueue.enqueue(readParams); 
            blockReadIsOver=0;  
        }
    }
    int leftMarker=overlapSNPcnt%BLOCK_SIZE;
    if(leftMarker>0){
        readParams.G=G.submat(0,0,leftMarker-1,gp_sampleNum-1);
        readParams.r=r.subvec(0,leftMarker-1);
        cout<<"["<< showtime() << "] ...Scanning " << lineCNT << " markers and overlapping markers " << overlapSNPcnt<< ".\n";
        foutLog<<"["<< showtime() << "] ...Scanning " << lineCNT << " markers and overlapping markers " << overlapSNPcnt<< ".\n";
        workQueue.enqueue(readParams); 
    }

    queueFilled = true;
    bool obtainedRes{false};
    int blockCnt=0;
    while ((obtainedRes = resultQueue.try_dequeue(kin)) or numWorking > 0) {
        if (obtainedRes) { break;}
    }
    for (auto& t:workers) { t.join(); }
    return kin;
}



void blockWorkerHomo(int & BLOCK_SIZE, int & gp_sampleNum, int &numWorking,
                    atomic<bool>& queueFilled,
                    moodycamel::ConcurrentQueue<readThreadParams>& workQueue,
                    moodycamel::ConcurrentQueue<fmat> & resultQueue){

    bool gotNewStart{false};
    double denom=0;

    fmat V1=ones <fmat> (gp_sampleNum, 1);
    readThreadParams m;
    fmat num=zeros <fmat> (gp_sampleNum, gp_sampleNum);

    while ((gotNewStart = workQueue.try_dequeue(m)) or !queueFilled) {
        if (!gotNewStart) {continue;}
        gotNewStart = false;
        m.G=0.5*m.G;
        fvec AF=mean(m.G,1);
        m.G -=AF*V1.t();
        int i=0;
        if(WEIGHT==1){
            fmat _num=m.G.t()*m.G;
            _num.diag()=sum(m.G%m.G%(m.r*V1.t()));
            num+=_num;
            denom+=sum(AF%(1-AF)%m.r%m.r);      
        }
        else if(WEIGHT==2){
            m.G/=(sqrt(AF%(1-AF))*V1.t());
            fmat _num=m.G.t()*m.G;
            _num.diag()=sum(m.G%m.G%(m.r*V1.t()));
            num+=_num;
            denom+=sum(m.r%m.r);
        }
    } 
    num/=denom;
    resultQueue.enqueue(num); 
    --numWorking;    
}

static fmat admixEst(string gp_file, string af_indiv_file, int gp_sampleNum, int BLOCK_SIZE, map<string, int> LIST, map <int, int> INDEX) {
    
    moodycamel::ConcurrentQueue <fmat> resultQueue;
    moodycamel::ConcurrentQueue <readThreadParams> workQueue;
    vector<thread> workers;
    atomic<bool> queueFilled{false};
    readThreadParams readParams;
    fmat kin=zeros <fmat> (gp_sampleNum, gp_sampleNum);
    fvec r = zeros <fvec> (BLOCK_SIZE);  
    fmat G = zeros <fmat> (BLOCK_SIZE, gp_sampleNum);    // dosage 
    fmat P = zeros <fmat> (BLOCK_SIZE, gp_sampleNum);    // individual AF 
    map<int, string> overlapSNPlist;
    int numWorking=1;
    workers.emplace_back(blockWorkerAdmix,ref(BLOCK_SIZE), ref(gp_sampleNum),ref(numWorking),ref(queueFilled),ref(workQueue),ref(resultQueue));

    VcfFileReader reader; 
    VcfHeader header;
    VcfRecord  record;
    reader.open(INPUT_FILE.c_str(),header);  

    VcfHeader af_indiv_header;
    VcfFileReader af_indiv_reader;
    VcfRecord  af_indiv_record;
    af_indiv_reader.open(AF_INDIV_FILE.c_str(),af_indiv_header);
    af_indiv_reader.readVcfIndex();

    int lineCNT=0;
    int blockReadIsOver=0;
    int matPtr=0;

    cout<<"["<< showtime() << "] Reading the markers in " << INPUT_FILE.c_str() << "\n";
    foutLog<<"["<< showtime() << "] Reading the markers in " << INPUT_FILE.c_str() << "\n";
    //omp_set_num_threads(50);
    while(reader.readRecord(record)){              
        VcfRecordInfo& info = record.getInfo();
        string chr = record.getChromStr();
        VcfRecordGenotype& geno = record.getGenotypeInfo();
        stringstream pos;
        pos << record.get1BasedPosition();
        string ID = chr+"_"+pos.str();

        if(LIST.count(ID)){
            double r2_tmp=1;
            double af_tmp=0.5;      
            fvec _G = zeros <fvec> (gp_sampleNum);  
            if(GENE_MODE.compare("DS")==0){
                r2_tmp = atof((*(info.getString("DR2"))).c_str());
                af_tmp = atof((*(info.getString("AF"))).c_str());
                int k=0;
                for (k=0; k<gp_sampleNum; k++){                   
                     _G[k]= atof((*geno.getString("DS",k)).c_str());
                }
            }
            else{ 
                af_tmp=0;    
                for (int k=0; k<gp_sampleNum; k++){   
                    string tmp= (*geno.getString("GT",k)).c_str();
                    _G[k]=genoCode[tmp];
                    af_tmp=af_tmp+genoCode[tmp]/2;
                }
                r2_tmp=1;
                af_tmp=af_tmp/gp_sampleNum;
            }
            
            if(r2_tmp >= R2_FILTER & af_tmp >= MAF_FILTER & af_tmp <= (1-MAF_FILTER)) {  //***
                overlapSNPcnt++;
                matPtr=overlapSNPcnt%BLOCK_SIZE;
                if(matPtr==0){
                    matPtr=BLOCK_SIZE;
                    blockReadIsOver=1;
                }

                overlapSNPlist[matPtr-1]=ID;
                G.row(matPtr-1)=_G.t();
                r[matPtr-1]=r2_tmp;

                vector<string> tokens = ::split(ID,string("_"));
                int pos_new = atof(tokens[1].c_str());
                af_indiv_reader.set1BasedReadSection(tokens[0].c_str(), pos_new-1, pos_new);
                af_indiv_reader.readRecord (af_indiv_record);
                VcfRecordGenotype& af_indiv = af_indiv_record.getGenotypeInfo();
                int k=0;
                //#pragma omp parallel for private(k)
                for (k=0; k<gp_sampleNum; k++){       
                    P(matPtr-1,k)= atof((*af_indiv.getString("AF1",INDEX[k])).c_str());
                    if(P(matPtr-1,k)<DELTA){P(matPtr-1,k)=DELTA;}
                    else if (P(matPtr-1,k)>1-DELTA){P(matPtr-1,k)=1-DELTA;}
                }
                if(blockReadIsOver){
                    cout<<"["<< showtime() << "] ...Scanning " << lineCNT << " markers and overlapping markers " << overlapSNPcnt<< ".\n";
                    foutLog<<"["<< showtime() << "] ...Scanning " << lineCNT << " markers and overlapping markers " << overlapSNPcnt<< ".\n";
                    readParams.G=G;
                    readParams.P=P;
                    readParams.r=r;
                    workQueue.enqueue(readParams); 
                    blockReadIsOver=0;  
                }
            }
        }
        lineCNT++;
    }

    int leftMarker=overlapSNPcnt%BLOCK_SIZE;
    if(leftMarker>0){
        readParams.G=G.submat(0,0,leftMarker-1,gp_sampleNum-1);
        readParams.P=P.submat(0,0,leftMarker-1,gp_sampleNum-1);
        readParams.r=r.subvec(0,leftMarker-1);
        cout<<"["<< showtime() << "] ...Scanning " << lineCNT << " markers and overlapping markers " << overlapSNPcnt<< ".\n";
        foutLog<<"["<< showtime() << "] ...Scanning " << lineCNT << " markers and overlapping markers " << overlapSNPcnt<< ".\n";
        workQueue.enqueue(readParams); 
    }

    queueFilled = true;
    bool obtainedRes{false};
    int blockCnt=0;
    while ((obtainedRes = resultQueue.try_dequeue(kin)) or numWorking > 0) {
        if (obtainedRes) { break;}
    }
    for (auto& t:workers) { t.join(); }
    return kin;
}




void blockWorkerAdmix(int & BLOCK_SIZE, int & gp_sampleNum, int &numWorking,
                    atomic<bool>& queueFilled,
                    moodycamel::ConcurrentQueue<readThreadParams>& workQueue,
                    moodycamel::ConcurrentQueue<fmat> & resultQueue){

    bool gotNewStart{false};
    fmat V1= ones <fmat> (gp_sampleNum, 1);
    readThreadParams m;
    fmat num=zeros <fmat> (gp_sampleNum, gp_sampleNum);
    fmat denom=zeros <fmat> (gp_sampleNum, gp_sampleNum);

    while ((gotNewStart = workQueue.try_dequeue(m)) or !queueFilled) {
        if (!gotNewStart) { continue; }
        gotNewStart = false; 
        m.G=0.5*m.G;
        m.G-=m.P%(m.r*V1.t());
        m.G-=mean(m.G,1)*V1.t();

        int i=0;
        if(WEIGHT==1){
            fmat adjustG=m.G+m.P%((m.r-sqrt(m.r))*V1.t());
            adjustG-=mean(adjustG,1)*V1.t();
            fmat _num=m.G.t()*m.G;
            _num.diag()=sum(adjustG%adjustG%(m.r*V1.t()));  // Adjust the self-kinship
            num+=_num;
            fmat H=sqrt(m.P%(1-m.P))%(m.r*V1.t());
            denom+=H.t()*H;
        }
        else if(WEIGHT==2){
            fmat adjustG=m.G+m.P%((m.r-sqrt(m.r))*V1.t());
            adjustG-=mean(adjustG,1)*V1.t();
            adjustG/=sqrt(m.P%(1-m.P));
            m.G/=sqrt(m.P%(1-m.P));
            fmat _num=m.G.t()*m.G;
            _num.diag()=sum(adjustG%adjustG%(m.r*V1.t()));  // Adjust the self-kinship
            num+=_num;
            fmat H=sqrt(m.r*V1.t());
            denom+=H.t()*H;
        }
    } 
    num/=denom;
    resultQueue.enqueue(num); 
    --numWorking; 
}


static map<string,int> getSNPlist(string inputfile){
    VcfRecord  record;
    VcfHeader  header;
    VcfFileReader reader;
    reader.open(inputfile.c_str(),header);
    map <string, int> LIST;
    int lineCNT=0;
    cout<<"["<< showtime() << "] Reading the markers in " << inputfile.c_str()<< "\n";
    foutLog<<"["<< showtime() << "] Reading the markers in " << inputfile.c_str()<< "\n";
    while(reader.readRecord(record)){                    
        string chr = record.getChromStr();
        stringstream pos;
        pos << record.get1BasedPosition();
        string ID = chr+"_"+pos.str();
        lineCNT++;
        LIST[ID]=lineCNT;
    }
    cout<<"["<< showtime() << "] Total " << lineCNT <<  " markers in " << inputfile.c_str()<< "\n";
    foutLog<<"["<< showtime() << "] Total " << lineCNT <<  " markers in " << inputfile.c_str()<< "\n";
    reader.close();
    return LIST;
}

static void output(fmat &KIN, string OUTPUT_FILE, int gp_sampleNum, string inputfile, int snp_count){

   
    ofstream FOUT;
    ofstream FOUT1;
    ofstream FOUT2;
    string tmp=OUTPUT_FILE+".kin";
    FOUT.open(tmp.c_str());
    tmp=OUTPUT_FILE+".ID";
    FOUT1.open(tmp.c_str());
    tmp=OUTPUT_FILE+".inbreed";
    FOUT2.open(tmp.c_str());
    cout<<"["<< showtime() << "] Write results to \"" << OUTPUT_FILE << ".*\"\n";
    FOUT << setiosflags(ios::fixed)<<setprecision(4);
    FOUT << "Ind1\tInd2\tNSNP\tkinship\n";
    FOUT2 << "IID\tINBREEDCOEF\n";       

    VcfHeader gp_header;
    VcfFileReader gp_reader;
    gp_reader.open(inputfile.c_str(),gp_header);
    FOUT2 << setiosflags(ios::fixed)<<setprecision(4);
    for(int i=0; i < gp_sampleNum; i++ ){
        FOUT1 << gp_header.getSampleName(i) <<"\n";
        FOUT2 << gp_header.getSampleName(i) <<"\t"<< 2*KIN(i,i)-1 <<"\n";
        if(i<gp_sampleNum-1){
            for(int j=i+1; j < gp_sampleNum; j++ ){
                FOUT << gp_header.getSampleName(i) << "\t" <<gp_header.getSampleName(j) << "\t" <<snp_count <<"\t"<< KIN(i,j)<< "\n";
            }  
        }    
    }
    FOUT.close();
    FOUT1.close();
    FOUT2.close();
    gp_reader.close();

    tmp=OUTPUT_FILE+".matrix";
    FOUT.open(tmp.c_str());
    FOUT << setiosflags(ios::fixed)<<setprecision(4);
    FOUT << KIN ;
    FOUT.close();
}

static int getVCFsampleNum(string inputfile){
    VcfHeader header;
    VcfFileReader reader;
    int sampleNum=0;
    if(reader.open(inputfile.c_str(),header)){
       sampleNum=header.getNumSamples();
       reader.close();
       return sampleNum;
    } 
    else{
        cout<<"["<< showtime() << "Error! Could not get the sample number of the genotype probability vcf file\n";
        exit(-1);
    }
}

static void checkHeader(int af_indiv_sampleNum, int gp_sampleNum, map<int,int> & INDEX, string gp_file, string af_indiv_file){
    map<string, int> af_indiv_sample_index;
    VcfHeader af_indiv_header;
    VcfFileReader af_indiv_reader;
    VcfHeader gp_header;
    VcfFileReader gp_reader;
    af_indiv_reader.open(af_indiv_file.c_str(),af_indiv_header);
    gp_reader.open(gp_file.c_str(),gp_header);
    af_indiv_reader.close();  
    gp_reader.close();
    for(int i=0; i<af_indiv_sampleNum; i++){
            af_indiv_sample_index[af_indiv_header.getSampleName(i)]=i;
    }
        
    for(int i=0; i<gp_sampleNum; i++){
        string tmp=gp_header.getSampleName(i);
        if(af_indiv_sample_index.count(tmp) !=0) {
            INDEX[i]=af_indiv_sample_index[tmp];
        }
        else{
            cout<<"["<< showtime() << "] Check the headers of individual-specific allele frequency file and the VCF files" << "\n";
            cout<<"["<< showtime() << "] the sample " << tmp << "  is not detected in the file" << af_indiv_file.c_str()<< "\n";
            exit (-1); 
        }
    }  
}

bool  display_usage(){
    fprintf ( stderr, "\nBasic usage: seekin kinship <command> [option] \n" );
    fprintf ( stderr, "  -i <str>     the input VCF file of studied samples (required)\n" );
    fprintf ( stderr, "  -o <str>     the prefix of output files \n" );
    fprintf ( stderr, "  -p <str>     the population structure mode: homo or admix (required)\n" );
    fprintf ( stderr, "  -w <int>     weighted function used to average all markers, 1: 2pqr^4; 2 2r^4 [1]\n" );
    fprintf ( stderr, "  -d <str>     use genotype or dosage GT or DS, [DS] \n" );
    fprintf ( stderr, "  -a <str>     the individual specific allele frequency file (required for admixture)\n");
    fprintf ( stderr, "  -r <float>   imputation dosage r2 filter [0.5] \n" );
    fprintf ( stderr, "  -m <float>   minor allele frequency filter [0.05] \n" );
    fprintf ( stderr, "  -n <int>     the number of SNPs per block [10,000]\n" );
    fprintf ( stderr, "  -t <int>     the number of threads [1]\n" );
    fprintf ( stderr, "  -h <>        help information\n\n" );
    return (true);
}

vector<string> split(const string& src, string sp) { 
    vector<string> strs; 
    int sp_len = sp.size(); 
    int position = 0, index = -1; 
    while (-1 != (index = src.find(sp, position))) { 
        strs.push_back(src.substr(position, index - position)); 
        position = index + sp_len; 
    } 
    string lastStr = src.substr(position); 
    if (!lastStr.empty()) 
        strs.push_back(lastStr); 
    return strs; 
} 

static void paraCheck (){

    int paraCheckFlag=1;
    int fileCheckFlag=1;

    ifstream fin;
    stringstream ss;
    ss << getpid();
    if (OUTPUT_FILE.compare(default_str)!=0) {
        LOG_FILE=OUTPUT_FILE;
    }

    genoCode["0|0"]=0;    // Need to change the defintion of this map !!!!!
    genoCode["0|1"]=1;
    genoCode["1|0"]=1;
    genoCode["1|1"]=2;
    genoCode["0/0"]=0;
    genoCode["0/1"]=1;
    genoCode["1/1"]=2;
    genoCode["./."]=0;

    LOG_FILE.append(".");
    LOG_FILE.append(ss.str());
    LOG_FILE.append(".log");
    foutLog.open(LOG_FILE.c_str());
    if(foutLog.fail()){
        cerr << "Error: cannot create the log file.\n";
    }


    cout << "********************\n";
    cout << "Module kinship\n";
    cout << "********************\n";
    foutLog << "********************\n";
    foutLog << "Module kinship\n";
    foutLog << "********************\n";
    fprintf ( stdout, "Parameters: \n" );
    fprintf ( stdout, " -i %s \n", INPUT_FILE.c_str() );
    fprintf ( stdout, " -a %s \n", AF_INDIV_FILE.c_str() );
    fprintf ( stdout, " -o %s \n", OUTPUT_FILE.c_str() );
    fprintf ( stdout, " -r %f \n", R2_FILTER);
    fprintf ( stdout, " -m %f \n", MAF_FILTER);
    fprintf ( stdout, " -d %s \n", GENE_MODE.c_str() );
    fprintf ( stdout, " -p %s \n", POP_STRUCT.c_str() );
    // The output of weight is always zero, fix me. 
    fprintf ( stdout, " -w %d \n", WEIGHT);
    fprintf ( stdout, " -n %d \n", BLOCK_SIZE);
    fprintf ( stdout, " -t %d \n\n", THREAD_NUM);


    //t1=omp_get_wtime();
    cout<<"["<< showtime() << "] Started!\n";
    foutLog<<"["<< showtime() << "] Started!\n";

    if (BLOCK_SIZE<1 || BLOCK_SIZE>1000000){
        cout<<"["<< showtime() << "] Error! The set block size ranges from 1 to 1000,000!\n";
        foutLog<<"["<< showtime() << "] Error! The set block size ranges from 1 to 1000,000!\n";
        display_usage();
        exit(-1);    
    }

    if (THREAD_NUM<1 || THREAD_NUM>100){
        cout<<"["<< showtime() << "] Error! The set thread ranges from 1 to 100!\n";
        foutLog<<"["<< showtime() << "] Error! The set thread ranges from 1 to 100!\n";
        display_usage();
        exit(-1);    
    }

    if (MAF_FILTER<0 || MAF_FILTER>=0.5){
        cout<<"["<< showtime() << "] Error! The MAF filter ranges from 0 to 0.5!\n";
        foutLog<<"["<< showtime() << "] Error! The MAF filter ranges from 0 to 0.5!\n";
        display_usage();
        exit(-1);    
    }

    if (R2_FILTER<0 || R2_FILTER>1){
        cout<<"["<< showtime() << "] Error! The R2 filter ranges from 0 to 0.5!\n";
        foutLog<<"["<< showtime() << "] Error! The R2 filter ranges from 0 to 0.5!\n";
        display_usage();
        exit(-1);    
    }



    if (POP_STRUCT.compare(default_str)==0){
        cout<<"["<< showtime() << "] Error! No population structure mode is specified!\n";
        foutLog<<"["<< showtime() << "] Error! No population structure mode is specified!\n";
        display_usage();
        exit(-1);    
    }
    else if(POP_STRUCT.compare("homo")!=0 & POP_STRUCT.compare("admix")!=0 ){
        cout<<"["<< showtime() << "] Error! The specified population structure mode is neither homo nor admix!\n";
        foutLog<<"["<< showtime() << "] Error! The specified population structure mode is neither homo nor admix!\n";
        display_usage();
        exit(-1);   
    }

    if(INPUT_FILE.compare(default_str)!=0){

        fin.open(INPUT_FILE.c_str());
        if(fin.fail()){
            cout << "["<< showtime() << "] Error! Fail to open the VCF file " << INPUT_FILE.c_str()  << "\n";
            foutLog << "["<< showtime() << "] Error! Fail to open the VCF file " << INPUT_FILE.c_str()  << "\n";
            exit(-1);
        }
        else{
            fin.close();
        }
    }
    else{
        cout << "["<< showtime() << "] No VCF file is specified, please check it!" << "\n";
        foutLog << "["<< showtime() << "] No VCF file is specified, please check it!" << "\n";
        exit(-1);
    }

    if(AF_INDIV_FILE.compare(default_str)!=0){
        fin.open(AF_INDIV_FILE.c_str());
        if(fin.fail()){
            cout << "["<< showtime() << "] Error! Fail to open the individual-specific allele frequency file " << AF_INDIV_FILE.c_str()  << "\n";
            foutLog << "["<< showtime() << "] Error! Fail to open the individual-specific allele frequency file " << AF_INDIV_FILE.c_str()  << "\n";

            exit(-1);
        }
        else{
            fin.close();
            indivAFlag=1;
        }
    }
    else{
        cout << "["<< showtime() << "] No individual-specific allele frequency file is specified!" << "\n";
        foutLog << "["<< showtime() << "] No individual-specific allele frequency file is specified!" << "\n";
        indivAFlag=0;
    }

    // determine the estimator mode 

    if(POP_STRUCT.compare("homo")==0 &  WEIGHT==2 ){
        cout << "["<< showtime() << "] Use the homo mode with weight function r^4" << "\n";
        foutLog << "["<< showtime() << "] Use the homo mode with weight function r^4" << "\n";
    }
    else if(POP_STRUCT.compare("homo")==0 & WEIGHT==1){
        cout << "["<< showtime() << "] Use the homo mode with weight function 2pqr^4" << "\n";
        foutLog << "["<< showtime() << "] Use the homo mode with weight function 2pqr^4" << "\n";
    }
    else if(POP_STRUCT.compare("admix")==0 & WEIGHT==2 & indivAFlag==1){
        cout << "["<< showtime() << "] Use the admix mode with weight function r^4" << "\n";
        foutLog << "["<< showtime() << "] Use the admix mode with weight function r^4" << "\n";
    }
    else if(POP_STRUCT.compare("admix")==0 & WEIGHT==1 & indivAFlag==1){
        cout << "["<< showtime() << "] Use the admix mode with weight function 2pqr^4" << "\n";
        foutLog << "["<< showtime() << "] Use the admix mode with weight function 2pqr^4" << "\n";
    }
    else{
        cout << "["<< showtime() << "] Error! No mode is identified, please check the input files" << "\n";
        foutLog << "["<< showtime() << "] Error! No mode is identified, please check the input files" << "\n";
        display_usage();
        exit(-1);
    }
    if(GENE_MODE.compare("DS")==0 |  GENE_MODE.compare("GT")==0 ){;}
    else{
        cout << "["<< showtime() << "] Error! No mode is identified, please use DS or GT" << "\n";
        foutLog << "["<< showtime() << "] Error! No mode is identified, please use DS or GT" << "\n";
        display_usage();
        exit(-1);
    }

    // determine the thread number used for kinship estimation 
    if(BLOCK_SIZE==default_int || BLOCK_SIZE==0){
        BLOCK_SIZE=10000;
        cout << "["<< showtime() << "] BLOCK_SIZE is set to the default value 1,000" << "\n";
        foutLog << "["<< showtime() << "] BLOCK_SIZE is set to the default value 1,000" << "\n";
    }
    else{
        cout << "["<< showtime() << "] BLOCK_SIZE is set to " << BLOCK_SIZE << "\n";
        foutLog << "["<< showtime() << "] BLOCK_SIZE is set to " << BLOCK_SIZE << "\n";
    }

    if(THREAD_NUM==default_int){
        THREAD_NUM=1;
        cout << "["<< showtime() << "] THREAD_NUM is set to the default value 1" << "\n";
        foutLog << "["<< showtime() << "] THREAD_NUM is set to the default value 1" << "\n";
    }
    else{
        cout << "["<< showtime() << "] THREAD_NUM is set to " << THREAD_NUM << "\n";
        foutLog << "["<< showtime() << "] THREAD_NUM is set to " << THREAD_NUM << "\n";
    }

}

static void initenv (int argc, char ** argv){
    char copt;
    int input=0;

    while((copt=getopt(argc, argv, "i:o:w:a:r:m:d:n:p:t:h")) != EOF){
        switch(copt){   
            case 'i':
                input = 1;
                INPUT_FILE=strdup(optarg);
                continue;
            case 'o':
                OUTPUT_FILE=strdup(optarg);
                continue;
            case 'w':
                WEIGHT=atof(optarg);
                continue;
            case 'd':
                GENE_MODE=strdup(optarg);
                continue;
            case 'a':
                AF_INDIV_FILE=strdup(optarg);
                continue;
            case 'r':
                R2_FILTER=atof(optarg);
                continue;
            case 'm':
                MAF_FILTER=atof(optarg);
                continue;
            case 'n':
                BLOCK_SIZE=atof(optarg);
                continue;
            case 'p':
                POP_STRUCT=strdup(optarg);
                continue;
            case 't':
                THREAD_NUM=atof(optarg);
                continue;
            case 'h':
                display_usage();
                exit(0);
        }
    } 
    if(input==0){
        display_usage();
        exit(-1);
    }
}

static string showtime(){
    time ( &rawtime );
    timeinfo = localtime ( &rawtime);
    string temp =asctime (timeinfo);
    return temp.substr(0,temp.size()-1);
}