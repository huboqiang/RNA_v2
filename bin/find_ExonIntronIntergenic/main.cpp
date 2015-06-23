// g++ -o main main.cpp -lboost_regex -lboost_iostreams

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

using namespace std;

void usage()
{
   cout << "./find_Exon_Intron_Intergenic <in.bed> <hg19.fa.fai> " << endl;
   cout << "  -h           get help information"   << endl;
   cout << "example: " << endl;
   cout << "./find_Exon_Intron_Intergenic /data/Analysis/huboqiang/project/human_embryo_sequencing/07.CHIP-seq/05.region-counts/region_definition/refseqGene.gtf.exon.sort.bed /data/Analysis/huboqiang/database/hg19/hg19.fa.fai " << endl;
   exit (0);
}

void load_str(string fai_file,map<string,string> &ChrStr){
	ifstream infile;
	infile.open(fai_file.c_str());
	if ( ! infile ){
		cerr << "fail to open input file" << fai_file << endl;
		exit(0);
	}
	string lineStr;
	while (getline(infile,lineStr,'\n')){
		if (lineStr[0] == ' ' || lineStr[0] == '\n'){
			continue;
		}
		vector<string> lineVec;
		boost::split(lineVec,lineStr,boost::is_any_of(":, \t\n"), boost::token_compress_on);
		int length =  boost::lexical_cast<int>(lineVec[1]);
		ChrStr[lineVec[0]].assign(length,'0');
	}
	infile.close();
}

void load_exon(string sam_file,map<string,string> &ChrStr){
	ifstream infile;
	string file_name = sam_file;
	infile.open(file_name.c_str());
	if ( ! infile ){
		cerr << "fail to open input file " << file_name << endl;
		exit(0);
	}
	string lineStr;
	
	// 1 for intron, 2 for Exon, 0 for intergenic
	string pre_chrom = "";
	while (getline(infile,lineStr,'\n')){
		if (lineStr[0] == ' ' || lineStr[0] == '\n'){
			continue;
		}
		vector<string> lineVec;
		boost::split(lineVec,lineStr,boost::is_any_of("\t\n"), boost::token_compress_on);
		int beg =  boost::lexical_cast<int>(lineVec[1]);
		int end =  boost::lexical_cast<int>(lineVec[2]);
		
		string chrom = lineVec[0];
		if (chrom != pre_chrom){
			cerr << "processing " << chrom << endl;
		}
		for ( int i=beg;i<end;i++){
			if ( ChrStr[chrom][i-1] == '0' ){
				ChrStr[chrom][i-1] = '1';	
			}
		}
		vector<string> l_beg;
		vector<string> l_end;
		boost::split(l_beg,lineVec[5],boost::is_any_of(","), boost::token_compress_on);
		boost::split(l_end,lineVec[6],boost::is_any_of(","), boost::token_compress_on);
		for ( int i=0;i<l_beg.size();i++ ){
			int exon_beg =  boost::lexical_cast<int>(l_beg[i]) + 1;
			int exon_end =  boost::lexical_cast<int>(l_end[i])    ;
			ChrStr[chrom].replace(exon_beg-1,exon_end-exon_beg+1,exon_end-exon_beg+1,'2');
		}
		pre_chrom = chrom;
	}
	infile.close();
}

int main(int argc, char *argv[])
{
   int c;
   while ( (c=getopt(argc,argv,"h")) != -1 ){
      switch(c)
      {
         case 'h' : usage();break;
         default : usage();
      }
   }
   if (argc < 2) usage();
	
	string sam_file   = argv[optind++];
	string fai_file   = argv[optind++];
	
	map<string,string> ChrStr;
	
	load_str(fai_file,ChrStr);
	load_exon(sam_file,ChrStr);
	
	map<string,string >::iterator  idx_chrom;
	string pre_char = "";
	string cnt_char = "";
	string stat     = "Intergenic";
		
	for ( idx_chrom=ChrStr.begin();idx_chrom!=ChrStr.end();idx_chrom++ ){
		string chrom = idx_chrom->first;
		int pos_beg = 0;
		int pos_end = 0;
		for ( int i=0;i<ChrStr[chrom].size();i++ ){
			cnt_char = ChrStr[chrom][i];
			if ( cnt_char != pre_char ){
			
				if (pre_char != ""){	
					cout << chrom << "\t" << pos_beg << "\t" << pos_end << "\t" << stat << endl;
					
					if (cnt_char == "0"){
						stat = "Intergenic";
					}
					if (cnt_char == "1"){
						stat = "Intron";
					}
					if (cnt_char == "2"){
						stat = "Exon";
					}

				}
				pos_beg = i+1;
				pos_end = i+1;
			}
			if ( cnt_char == pre_char ){
				pos_end = i+1;
			}
			pre_char = cnt_char;	
		}
		cout << chrom << "\t" << pos_beg << "\t" << pos_end << "\t" << stat << endl;
	} 
}