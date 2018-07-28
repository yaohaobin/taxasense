#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<stack>
#include<algorithm>
#include<sstream>
#include<stdlib.h>

#include <sdsl/suffix_trees.hpp>

using namespace std;
using namespace sdsl;




void extract_diff(string& large, string& small,string outdir){
    cst_sct3<csa_bitcompressed<>> cst;

    string query = large + "$" + small;

    int_vector<2> belong(cst.size(),0);

    for(uint32_t i = large.length();i<cst.csa.size();i++)
    	belong[i]=1;

    construct_im(cst,query,1);

    int_vector<8>table(cst.csa.size(),0);
    for(unsigned long int i = 1;i<cst.csa.size();i++)
             table[i] = belong[cst.csa[i]];

    belong.resize(0);
    wt_huff<>wt;
    construct_im(wt,table);
    table.resize(0);
    cout<<"wt sigma: "<<wt.sigma<<" wt size: "<<wt.size()<<endl;
    double mean = 0;
    int num = 0;
    int max = 0;
    uint32_t lastnodedepth = 0,lastpos = 0,lastdepth=0;
     
    int_vector<8>common(cst.csa.size(),0);

    for (auto it=cst.begin(); it!=cst.end(); ++it) {
    	auto v = *it;
    	if(cst.node_depth(v)<=lastnodedepth){
                if(lastdepth != 0){
                   if(common[lastpos] < lastdepth) common[lastpos] = lastdepth;
                }
                lastdepth = lastpos = 0;
                
            }
        lastnodedepth = cst.node_depth(v);
        uint8_t idx = 1;
        uint64_t rankr = wt.rank(cst.rb(v),idx);
        uint64_t rankl = wt.rank(cst.lb(v),idx);



        if(rankr == rankl){
            lastdepth = cst.depth(v);
            lastpos = cst.csa[cst.lb(v)] + lastdepth - 1;
        }

        if(cst.depth(v) >=100 )
        	it.skip_subtree();


    }
    if(lastdepth != 0){
        if(common[lastpos] < lastdepth) common[lastpos] = lastdepth;
    }
    string outfile = outdir+".uni.fa";
    ofstream fout(outfile.c_str());
    uint64_t startpos;
    for(uint64_t i = 0;i<common.size();i++)
        if(common[i] != 0){
            startpos = i + 1 - common[i];
         
            fout<<seqtext.substr(startpos,common[i])+"$"<<endl;
            if(common[i] > max)
                max = common[i];
            mean = (mean*num+common[i])/(num+1);
            num++;
       }
    fout.close();
    cout<<"mean: "<<mean<<" num "<<num<<" max "<<max<<endl;



    




}

void extractdiff(string parentdir,vector<string>& dirs,vector<string>& childrentax){
    ifstream parentfile(parentdir.c_str());
    string nextline;
    string parentseq = "";
    getline(parentfile,nextline);
    while(parentfile>>nextline)
    	parentseq+=nextline;
    parentfile.close();
    cout<<"parentseq len: "<<parentseq.length()<<endl;

    for(unsigned i = 0;i<dirs.size();i++){
        string filename = dirs[i];
        ifstream seqfile(filename.c_str());
        string line;
        getline(seqfile,line);
        string seq = "";
        while(seqfile>>line)
            seq+=line;
        seqfile.close();
        cout<<"childseq len: "<<seq.length()<<endl;    
        extract_diff(seq,parentseq,childrentax[i]);    
        
     }
        
}