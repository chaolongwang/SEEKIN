/*
 *  Copyright (C) 2015-2016  Genome institute of Singapore (GIS),
 *                           Jinzhuang Dou and Chaolong Wang
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
using namespace std;
using namespace arma;
#include <boost/algorithm/string.hpp>

//*************************************************


const string default_str = "---this-is-a-default-string---";
const int default_int = 20;
const double default_double = -9.99999999;

string inputStudyCoord = default_str;
string inputRefBeta = default_str;
string outputStudyAF =default_str;
static string outputLog="SEEKIN_pcRelateAF.log";


int pcaCnt=default_int;
static time_t t1,t2,rawtime;
static float runningtime;
static struct tm * timeinfo;


bool getAF_display_usage();


static void initenv ( int argc, char ** argv );
static void paraCheck ();
static void getPCrelateAF ();
static vector <string> split(const string& src, string sp);
static string showtime();



int getAF(int argc, char ** argv){

    // check the input parameters and set default values. 
    initenv (argc, argv);
    // check the run mode
    paraCheck();
    // estimate the PC-related coefficients 
    getPCrelateAF();
}






static void getPCrelateAF(){



	//Reading the PCA coordinatefile of the studied samples

	map <int, string> studyId;

	cout<<"["<< showtime() << "] Read the PC coordinate file of the studied samples... " << "\n";
		
		
	int studyPCnum=0;
	int cnt;
	string line;
	int flag=0;
	ifstream FIN(inputStudyCoord.c_str());
	if(FIN.is_open()){
			cnt=0;
			getline(FIN,line);
			vector<string> tokens = ::split(line.c_str(),string("\t"));
			flag=0;
			for(int i=0;i<tokens.size();i++){
				if(tokens[i].compare("PC1")==0){
					flag=i;
					cout<<"["<< showtime() << "] the PC1 is in the " << flag << "th column [0-based]" << "\n";  
				}
    		}
    		if(flag==0){
    			cout<<"["<< showtime() << "] No PC1 flag  is detected, please check the format of the PC coordinate file ! << \n";
    			exit(-1);  
    		}
			studyPCnum= tokens.size()-flag;
			while(getline(FIN,line)){cnt++;}
	}
	FIN.close();

	int studyIDnum=cnt;
	mat studyPCoord = zeros <mat>(cnt,pcaCnt+1);

    cnt=0;
    FIN.open(inputStudyCoord.c_str());
	if(FIN.is_open()){
			getline(FIN,line);
			vector<string> tokens = ::split(line.c_str(),string("\t"));   // One debug about the split string here! 
			// Get the col of PC1 flag
			if(tokens.size()-flag> pcaCnt){
				cout<<"["<< showtime() << "] The PCs[" << tokens.size()-flag <<"]in studied coordinate file is larger than the specified PC numbers [" << pcaCnt <<"]\n"; 
				
			}
			else if(tokens.size()-flag < pcaCnt) {
				cout<<"["<< showtime() << "] Error: the PCs[" << tokens.size()-flag <<"]in studied coordinate file is smaller than the specified PC numbers [" << pcaCnt <<"]\n"; 
				exit(-1);
			}

			std::string delimiters("\t");
			while(getline(FIN,line)){
				cnt++;
				tokens = ::split(line.c_str(),string("\t"));
				studyId[cnt-1]= tokens[1];
				studyPCoord(cnt-1,0)=1;	
				for(int i=flag; i<=flag-1+pcaCnt; i++){  
					studyPCoord(cnt-1,i-flag+1)=atof(tokens[i].c_str());
				}
			}
	}
	FIN.close();
    cout<<"["<< showtime() << "] Finish reading the PC coordinate file of the studied sample numbers: "<< studyIDnum << "\n"; 


    // estimate the individual-specific allele frequencies 

    string newName=outputStudyAF+".gz";
    IFILE OUTaf = ifopen(newName.c_str(),"wb",InputFile::BGZF);
    ifprintf(OUTaf, "##fileformat=VCFv4.2\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequencies averaged across all individuals\">\n");
    ifprintf(OUTaf, "##FORMAT=<ID=AF1,Number=A,Type=Float,Description=\"Estimated individual-specific allele frequencies\">\n");
    ifprintf(OUTaf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for(int i=0;i<studyIDnum;i++){
    	ifprintf(OUTaf, "\t%s", studyId[i].c_str());
    }
    ifprintf(OUTaf, "\n");



	FIN.open(inputRefBeta.c_str());
	cnt=0;
	if(FIN.is_open()){
		getline(FIN,line);

		while(getline(FIN,line)){
			cnt++;

			if(cnt% 5000 ==0){ 
				  cout<<"["<< showtime() << "] " << cnt << " marker'allele frequencies have been calculated...\n";
			}
			mat beta = zeros <mat>(pcaCnt+1,1);
			vec af = zeros <vec>(studyIDnum);
			vector<string> tokens = ::split(line.c_str(),string("\t"));
			for(int i=0;i<=pcaCnt;i++){beta(i)=atof(tokens[i+5].c_str());}
			af=0.5*studyPCoord*beta;
			ifprintf(OUTaf, "%s\t%s\t.\t%s\t%s\t.\t.\tAF=%.4f\tAF1",tokens[0].c_str(),tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str(),mean(af));

			for(int i=0;i<studyIDnum;i++){
				// Fix the boundary 
				if(af(i) < 0.001) {af(i) = 0.001;}
				if(af(i) > 0.999) {af(i) = 0.999;}
				ifprintf(OUTaf, "\t%.4f", af(i)); 
			}
        	ifprintf(OUTaf, "\n");
        }
	}
	FIN.close();
	ifclose(OUTaf);
	t2 = clock();
	runningtime = (t2-t1)/CLOCKS_PER_SEC;
	cout<<"["<< showtime() << "] Finished!\n";
	cout << "Total CPU time: " << runningtime << " seconds.\n" << endl;
}


bool getAF_display_usage ()
{
	cout << R"(
Usage: seekin getAF [option]

Options:
	-i  File name of coordinates of study individuals in the reference space. [no default]
	-b  File name of the allele frequency model. [default "AF.model"]
	-k  Number of PCs used to compute allele frequencies. [default 2]
	-o  Output file name of individual-specific allele frequencies. [default "indivAF.vcf.gz"]

Note: the coordinate file (-i) can be prepared using the LASER software 
      (http://csg.sph.umich.edu/chaolong/LASER/)
Example: seekin getAF  -i study.coord  -b ref.AF.model -k 4 -o indivAF

)";
  return (true);
}



static void paraCheck (){

	int paraCheckFlag=1;
	int fileCheckFlag=1;
	ifstream fin;
	
	if (outputStudyAF.compare(default_str)!=0) {
		outputLog=outputStudyAF+".log";
	}

	//freopen(outputLog.c_str(),"w",stdout);


	if (outputStudyAF.compare(default_str)==0) {
		cout << "\nError: No individual-specific allele frequencies file is specified! \n";
		paraCheckFlag=0;
	}

	if (inputStudyCoord.compare(default_str)==0) {
		cout << "\nError: No study coordinates file is specified! \n";
		paraCheckFlag=0;
	}

	if (inputRefBeta.compare(default_str)==0) {
		cout << "\nError: No PC-related coefficients file is specified! \n";
		paraCheckFlag=0;
	}


	if(paraCheckFlag==0){
		getAF_display_usage();
		exit(-1);
	}



	cout << "********************\n";
	cout << "getAF\n";
	cout << "********************\n";


	t1=clock();
	cout<<"["<< showtime() << "] Started\n";

	fprintf ( stdout, "Parameters: \n" );
	fprintf ( stdout, " -i %s \n", inputStudyCoord.c_str() );
	fprintf ( stdout, " -b %s \n", inputRefBeta.c_str() );
	fprintf ( stdout, " -o %s \n", outputStudyAF.c_str() );
	fprintf ( stdout, " -k %d \n\n",pcaCnt );


	// check the format of input files 

	fin.open(inputStudyCoord.c_str());
	if(fin.fail()){
		cout << "Error: cannot find the study coordinates file '" << inputStudyCoord.c_str() << "'.\n" ; 
		fileCheckFlag=0;     
	}
	fin.close();

	fin.open(inputRefBeta.c_str());
	if(fin.fail()){
		cout << "Error: cannot find the PC-related coefficients file  '" << inputRefBeta.c_str()  << "'.\n" ; 
		fileCheckFlag=0;     
	}
	fin.close();

	if(fileCheckFlag==0){
		getAF_display_usage();
		exit(-1);
	}
}



static void initenv (int argc, char ** argv){
	char copt;
	int paraNum=0;
	while((copt=getopt(argc, argv, "i:b:o:k:h")) != EOF){
		switch(copt){   
			case 'i':
				inputStudyCoord=strdup(optarg);
				paraNum++;
				continue;
			case 'b':
				inputRefBeta=strdup(optarg);
				paraNum++;
				continue;
			case 'o':
				outputStudyAF=strdup(optarg);
				paraNum++;
				continue;
			case 'k':
				pcaCnt=atof(optarg);
				paraNum++;
				continue;
			case 'h':
				getAF_display_usage();
				exit(0);
		}
	} 
	if(paraNum==0){
		getAF_display_usage();
		exit(-1);
	}
}

static vector <string> split(const string& src, string sp) { 
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


static string showtime(){
	time ( &rawtime );
	timeinfo = localtime ( &rawtime);
	string temp =asctime (timeinfo);
	return temp.substr(0,temp.size()-1);
}