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
#include <stdio.h>

using namespace std;
using namespace arma;


/*************************************************
Function:
	modelAF.cpp
Description:
	The funciton for estimating the regression coefficients for the specified reference panel. 
Usage:
    @see  modelAF_display_usage()

*************************************************/


const string default_str = "---this-is-a-default-string---";
const int default_int = 20;
const double default_double = -9.99999999;

string inputRefVcf = default_str;
string inputRefCoord = default_str;
string outputRefBeta=default_str;
static string outputLog="SEEKIN_modelAF.log";
int pcaNum=default_int;

static time_t t1,t2,rawtime;
static float runningtime;
static struct tm * timeinfo;





bool modelAF_display_usage();

static vector <string> split(const string& src, string sp);
static void initenv ( int argc, char ** argv );
static void paraCheck ();
static void getPCrelate ();
static string showtime();



int modelAF(int argc, char ** argv){



    initenv (argc, argv);
    paraCheck();
    getPCrelate(); 
  	fclose(stdout);


}

//*****************************************
//    estimate the PC-related coefficients 
//*****************************************
static void getPCrelate(){


	int refPCnum=0;
	int cnt=0;
	int flag=0;
	string line;
	map <string, int> refPCoordId;	
	cout<<"["<< showtime() << "] Started...\n";
	ifstream FIN(inputRefCoord.c_str());
	if(FIN.is_open()){
			getline(FIN,line);
			vector<string> tokens = ::split(line.c_str(),string("\t"));
			flag = 0;
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
			refPCnum= tokens.size()-flag;
			while(getline(FIN,line)){cnt++;}
	}
	FIN.close();
	int refIDnum=cnt;
	if(pcaNum>refPCnum){pcaNum=refPCnum;}
	mat refPCoord = zeros <mat>(cnt,pcaNum+1);
    cnt=0;
    FIN.open(inputRefCoord.c_str());
	if(FIN.is_open()){
			getline(FIN,line);
			while(getline(FIN,line)){
				cnt++;	
				vector<string> tokens = ::split(line.c_str(),string("\t"));
                refPCoordId[tokens[1].c_str()]=cnt-1;
				refPCoord(cnt-1,0)=1;
				for(int i=flag; i<= flag + pcaNum -1 ; i++){		
					refPCoord(cnt-1,i-flag+1)=atof(tokens[i].c_str());     // Note the coord file format from LASER !
				}
			}
	}
	FIN.close();

    cout<<"["<< showtime()<< "] " << refIDnum << " individuals and "<< refPCnum <<" PCs are detected in the reference coordinate file\n"; 


	ofstream OUTbeta;

	OUTbeta.open(outputRefBeta.c_str());
	//OUTbeta << setiosflags(ios::fixed)<<setprecision(4);
    OUTbeta << "CHROM\tPOS\tREF\tALT\tAF";
    for(int i=0;i<=pcaNum;i++){OUTbeta << "\tbeta"<<i;}
    OUTbeta <<"\n";


	// Open the vcf file for reading.

    VcfFileReader reader;
    VcfHeader header;
    VcfRecord record;
    ofstream OUTaf;
    reader.open(inputRefVcf.c_str(),header);
    
	int sampleNum=header.getNumSamples();
    cout<<"["<< showtime() << "] " << sampleNum << " individuals are detected in the reference vcf file\n";


	// Check the sample IDs between the ref vcf file and  refcoord file

	map<int,string> vcfSampleID;
	vec W = zeros <vec>(sampleNum);
	for(int i=0; i<sampleNum; i++){
        	vcfSampleID[i]=header.getSampleName(i);
        	if(refPCoordId.count(vcfSampleID[i]) != 0){
        		W[i]=refPCoordId[vcfSampleID[i]];
        	}
        	else{
        		W[i]=-9;
        		cout<<"["<< showtime() << "] " << vcfSampleID[i] << " is not detected in " << inputRefCoord.c_str()<< "\n";
        	}
        	
    }

	int vcfLine=0;
    while(reader.readRecord(record)){     
        	vcfLine++;
        	if(vcfLine% 5000 ==0){
        		cout<<"["<< showtime()<< "] " << vcfLine << " markers' PC-related coefficients have been calculated\n";
        	}
    		string chr = record.getChromStr();
			stringstream pos;
			pos << record.get1BasedPosition();
    		VcfRecordGenotype& geno = record.getGenotypeInfo();
    		VcfRecordInfo& info = record.getInfo();
    		double afAll= atof((*(info.getString("AF"))).c_str());
    		string altallele = record.getAltStr();
    		string refallele = record.getRefStr();
    		mat GT = zeros <mat>(refIDnum,1);
    		mat error = zeros <mat>(sampleNum+1,1);
    		mat beta = zeros <mat>(pcaNum+1,1);
    		int i;

			for (i=0; i<sampleNum; i++){
				if(W[i]>-1){
					const string str = *geno.getString("GT",i);
					int _tmp=str.at(0)-48;
					if(_tmp>=0 && _tmp<=2){
						GT(W[i]) = _tmp+str.at(2)-48; // ASCII char to int
					}
					else{GT(W[i])=-9;}
				}
			}

        	// The mean was assigned to the missing genotype
			uvec mis = find(GT==-9);
			if(mis.n_elem>0){ 
            	double Gm = mean(GT(find(GT!=-9)));
            	for(i=0; i<mis.n_elem; i++){GT(mis(i)) = Gm; }
        	}
        
        	// Estimate beta using the Linear least squares method.
        	beta=inv(refPCoord.t()*refPCoord)*refPCoord.t()*GT;

       		OUTbeta << chr <<"\t"<<pos.str()<<"\t"<< refallele <<"\t"<<altallele<<"\t" << afAll;
        	for(int i=0;i<=pcaNum;i++){OUTbeta << "\t"<<beta(i);}
			OUTbeta <<"\n";
    }
    OUTbeta.close();
    t2 = clock();
	runningtime = (t2-t1)/CLOCKS_PER_SEC;
	cout<<"["<< showtime() << "] Finished!\n";
	cout << "Total CPU time: " << runningtime << " seconds.\n" << endl;

}

//*****************************************
//    dsiplay_usage  
//*****************************************


bool modelAF_display_usage ()
{
	cout << R"(
Usage: seekin modelAF [option]

Options:
	-i  File name of genotypes of reference individuals (gzipped VCF). [no default]
	-c  File name of PCA coordinate of reference individuals. [no default]
	-k  Number of PCs used to model allele frequencies. [default 2]
	-o  Output file name. [default "AF.model"]

Note: the coordinate file (-c) can be prepared using the LASER software 
      (http://csg.sph.umich.edu/chaolong/LASER/)
Example: seekin modelAF  -i reference.vcf.gz  -c ref.coord -k 4 -o ref.AF.model

)";
  return (true);
}



//*****************************************
//    split the string  
//*****************************************
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


//*****************************************
//    showtime  
//*****************************************
static string showtime(){
	time ( &rawtime );
	timeinfo = localtime ( &rawtime);
	string temp =asctime (timeinfo);
	return temp.substr(0,temp.size()-1);
}


//*****************************************
//    paraCheck  
//*****************************************
static void paraCheck (){

	int paraCheckFlag=1;
	int fileCheckFlag=1;
	ifstream fin;

	if (outputRefBeta.compare(default_str)!=0) {
		outputLog=outputRefBeta+".log";
	}
	//freopen(outputLog.c_str(),"w",stdout);
	cout << "********************\n";
	cout << "Module modelAF\n";
	cout << "********************\n";

	t1=clock();
	
 
	if (outputRefBeta.compare(default_str)==0) {

		cout << "\nError: No output PC-related coefficients file is specified! \n";
		paraCheckFlag=0;
	}

	if (inputRefVcf.compare(default_str)==0) {
		cout << "\nError: No reference vcf file is specified! \n";
		paraCheckFlag=0;
	}

	if (inputRefCoord.compare(default_str)==0) {
		cout << "\nError: No reference coordinates file is specified! \n";
		paraCheckFlag=0;
	}



	if(paraCheckFlag==0){
		modelAF_display_usage();
		exit(-1);
	}

	// print the parameters
	
	fprintf ( stdout, "Parameters: \n" );
	fprintf ( stdout, " -i %s \n", inputRefVcf.c_str() );
	fprintf ( stdout, " -c %s \n", inputRefCoord.c_str() );
	fprintf ( stdout, " -o %s \n", outputRefBeta.c_str() );
	fprintf ( stdout, " -k %d \n\n", pcaNum );




	// check the format of input files 

	fin.open(inputRefVcf.c_str());
	if(fin.fail()){
		cout << "Error: cannot find the reference vcf file  '" << inputRefVcf.c_str()  << "'.\n" ; 
		fileCheckFlag=0;     
	}
	fin.close();

	fin.open(inputRefCoord.c_str());
	if(fin.fail()){
		cout << "Error: cannot find the reference coordinates file  '" << inputRefCoord.c_str()  << "'.\n" ; 
		fileCheckFlag=0;     
	}
	fin.close();

	if(fileCheckFlag==0){
		modelAF_display_usage();
		exit(-1);
	}
}



static void initenv (int argc, char ** argv){
	char copt;
	int paraNum=0;
	while((copt=getopt(argc, argv, "i:c:o:k:h")) != EOF){
		switch(copt){   
			case 'i':
				inputRefVcf=strdup(optarg);
				paraNum++;
				continue;
			case 'c':
				inputRefCoord=strdup(optarg);
				paraNum++;
				continue;
			case 'o':
				outputRefBeta=strdup(optarg);
				paraNum++;
				continue;
			case 'k':
				pcaNum=atof(optarg);
				paraNum++;
				continue;
			case 'h':
				modelAF_display_usage();
				exit(0);
		}
	} 
	if(paraNum==0){
		modelAF_display_usage();
		exit(-1);
	}
}

