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

int gentext(int_vector<8>& seqtext,string& seq){
    int num = 0;
    unsigned long insertpos = seqtext.size();
    seqtext.resize(seqtext.size()+seq.size()+1);


    uint8_t text=0,revcomp=0,temp=0;
    int idx = 0;


    for(size_t i=0;i<seq.length();i++){



        switch(seq[i]){
            //case 'a': text[i] = 1;break;
            case 'A': temp = 0;break;
            case 'C': temp = 1;break;
            case 'T': temp = 2;break;
            case 'G': temp = 3;break;
            default: temp = 4;num++;

        }

                if(temp<4){
                        text= (text<<2) | temp;
                        //revcomp = revcomp | temp<<(2*idx);

                        idx++;
                        if(idx>2){
                                seqtext[insertpos] =text + 1;
                                //seqtext[insertpos] =text<revcomp?text:revcomp + 1;
                                insertpos++;
                                text = 0;
                                revcomp = 0;
                                idx = 0;
                        }
                }



    }

        if(idx !=0){
                seqtext[insertpos] = text + 1;
            insertpos++;
        }
    seqtext[insertpos] = 65;
    if(seqtext.size()>(insertpos+1) )
        seqtext.resize(insertpos+1);


   return num;
}



void getresultstr(vector<string>& seqs,string outdir,string sp,int_vector<8>& resultstr){
    

    string seqtext = "";
    for(unsigned i = 0;i<seqs.size();i++){
                //if(seqs[i].length() <10000)continue;
                //gentext(seqtext,seq);
                seqtext+=seqs[i];
                
                //seqtext.resize(seqtext.size()+1);
                seqtext+='Z';
                //seqtext[seqtext.size()-1] = 66;

     }

     string commonfile = outdir+'/'+sp+".sp";
     ifstream common(commonfile.c_str());
     uint64_t startpos = 0,endpos = 0;
     unsigned int depth = 0;
     while(common>>endpos>>depth){
         startpos = endpos + 1 - depth;
         
         
         string str= seqtext.substr(startpos,depth);

         gentext(resultstr,str);          
         
     }
     common.close();
    
     cout<<sp<<" "<<resultstr.size()<<endl;

}







void gencommon(vector<string>& seqs,string outdir,string taxid){
    wt_huff<>wt;
    int_vector<8>belong;
    string seqtext="";
    for(unsigned i = 0;i<seqs.size();i++){
                //if(seqs[i].length() <10000)continue;
                //gentext(seqtext,seq);
                seqtext+=seqs[i];
                unsigned long oldsize= belong.size();
                belong.resize(seqtext.length()+1);
                for(unsigned long j = oldsize;j<belong.size();j++)
            belong[j] = i;
                //seqtext.resize(seqtext.size()+1);
                seqtext+='Z';
                //seqtext[seqtext.size()-1] = 66;

     }
     cout<<seqtext.length()<<" "<<seqs.size()<<endl;
     if(seqtext.length()<100000)return;
     cst_sct3<csa_bitcompressed<>> cst;
     construct_im(cst,seqtext,1);
 

     int_vector<8>table(cst.csa.size(),0);
     for(unsigned long int i = 1;i<cst.csa.size();i++)
             table[i] = belong[cst.csa[i]];

     belong.resize(0);
    
     construct_im(wt,table);
     cout<<"wt sigma: "<<wt.sigma<<" wt size: "<<wt.size()<<" table[1] "<<(int)table[10000]<<" wt[1] "<<(int)wt[10000]<<endl;

     

     double mean = 0;
     int num = 0;
     int max = 0;
     uint32_t lastnodedepth = 0,lastpos = 0,lastdepth=0;
     
     int_vector<8>common(cst.csa.size(),0);
     for (auto it=cst.begin(); it!=cst.end(); ++it) {

        if (it.visit() == 1) {  // node visited the first time
            auto v = *it;    
            if(cst.node_depth(v)<=lastnodedepth){
                if(lastdepth != 0){
                   if(common[lastpos] < lastdepth) common[lastpos] = lastdepth;
                }
                lastdepth = lastpos = 0;
                
            }
            lastnodedepth = cst.node_depth(v);
            if(cst.depth(v)<4)
                continue;
            else{
                
                
               // get the node by dereferencing the iterator
              // if depth node is <= max_depth
                // process node, e.g. output it in format d-[lb, rb]
               
                //cout<<cst.depth(v)<<" "<<cst.lb(v)<<" "<<cst.rb(v)<<" leaf num: "<<cst.size(v)<<" "<<cst.lb(cst.select_leaf((int)cst.size(v)/2+1))<<" "<<cst.rb(cst.select_leaf((int)cst.size(v)))<<endl;
                //cout<<cst.depth(v)<<"-["<<cst.lb(v)<< ","<<cst.rb(v)<<"]"<<endl;
                uint8_t idx = 0;
                
                for(idx=0;idx<wt.sigma;idx++){
                    uint64_t rankr = wt.rank(cst.rb(v),idx);
                    uint64_t rankl = wt.rank(cst.lb(v),idx);
                    if(rankr == rankl){
                        it.skip_subtree();
                        break;
                    }    

                }
                if(idx == wt.sigma){
                   
                   lastdepth = cst.depth(v);
                   lastpos = cst.csa[cst.lb(v)] + lastdepth - 1;
                   

                   if(cst.depth(v)>=100)
                       it.skip_subtree();
                }
            }
        }

    }
    if(lastdepth != 0){
        if(common[lastpos] < lastdepth) common[lastpos] = lastdepth;
    }
    string outfile = outdir+".fa";
    ofstream fout(outfile.c_str());
    fout<<">"<<taxid<<endl;

    uint64_t startpos = 0;
    for(uint64_t i = 0;i<common.size();i++)
        if(common[i] != 0){
            startpos = i + 1 - common[i];
         
            fout<<seqtext.substr(startpos,common[i])<<endl;
            if(common[i] > max)
                max = common[i];
            mean = (mean*num+common[i])/(num+1);
            num++;
       }
    fout.close();
    cout<<"mean: "<<mean<<" num "<<num<<" max "<<max<<endl;


    

}




void gencommon(vector<string>&seqs,string dep){
    cst_sct3<csa_bitcompressed<>> cst;
	string seqtext = "";
	int_vector<8>belong;
    for(unsigned i = 0;i<seqs.size();i++){
        
		//gentext(seqtext,seq);
		seqtext+=seqs[i];
		unsigned long oldsize= belong.size();
		belong.resize(seqtext.length()+1);
		for(unsigned long j = oldsize;j<belong.size();j++)
            belong[j] = i;		
		//seqtext.resize(seqtext.size()+1);
		seqtext+='Z';		
		//seqtext[seqtext.size()-1] = 66;		
	
	}    
	cout<<seqtext.length()<<" "<<seqs.size()<<endl;
	construct_im(cst,seqtext,1);
	cout<<size_in_mega_bytes(cst)<<" "<<cst.csa.sigma<<endl;
	vector<unsigned int>common(seqtext.length(),0);
	
	map<unsigned long,unsigned int>coverpos;
	set<unsigned long>coverseq;
	
	unsigned int cur_depth = 0;
	bool prune = false;
	
	vector<pair<unsigned long,unsigned long> > intervals;
	
	unsigned long lb=1;
	
	for(size_t i = 2;i<cst.lcp.size();i++){
             

             //unsigned long endpos = cst.csa[i-1]+ cst.lcp[i] -1;
             /*
             if(common[endpos] >= cst.lcp[i]){
                    lb = i;
                    cur_depth = 0;
                    continue;
             }
             else
                   common[endpos] = cst.lcp[i];
             */
	     if(cur_depth == 0){					 
		     cur_depth = cst.lcp[i];
			 if(cur_depth !=0)
			     lb = i-1;
			 else
			     lb = i;
		 }
		 else{
		     if(cur_depth != cst.lcp[i]){
			     if( i-lb>=seqs.size() )
			         intervals.push_back(make_pair(lb,i-1));
				 cur_depth = cst.lcp[i];
				 lb = i-1;
			 }
		}
                
               
			     
	}
	if(cur_depth>0 && cst.lcp.size()-lb>=seqs.size())
	    intervals.push_back(make_pair(lb,cst.lcp.size()-1));
	cout<<"interval: "<<intervals.size()<<endl;	
	vector<pair<unsigned long,unsigned long> >commonseed;
	for(size_t i=0;i<intervals.size();i++){
	    coverseq.clear();
		for(size_t j=intervals[i].first;j<=intervals[i].second;j++)
		    coverseq.insert(belong[cst.csa[j]]);
		if(coverseq.size() == seqs.size()){
	            commonseed.push_back(intervals[i]);
                    //cout<<cst.csa[intervals[i].first]<<" "<<cst.csa[intervals[i].second]<<endl;
                }
	}
	cout<<"commonseed: "<<commonseed.size()<<endl;
	intervals.clear();
	bool repeat = false;

        
        double meanlen = 0;
        int num=0;  	
        string resultstr="";
	for(size_t i=0;i<commonseed.size();i++){
	    prune = false;
		cur_depth = cst.lcp[commonseed[i].first + 1];
                
	    for(size_t j=commonseed[i].first;j<=commonseed[i].second;j++){
		    unsigned long endpos = cst.csa[j] + cur_depth - 1;
			if(cur_depth < common[endpos] ) 
			    prune = true;
			else{
			    common[endpos] = cur_depth;
			}
		}
		if(!prune){
		    intervals.push_back(commonseed[i]);
                    meanlen = (meanlen*num+cur_depth)/(num+1);
                    num++;
                    //cout<<commonseed[i].first<<" "<<commonseed[i].second<<" "<<cur_depth<<endl;
                }
	}


		
	cout<<"final inter: "<<intervals.size()<<" "<<meanlen<<endl;	
    	
	
      


}

void extract(vector<string>&dirs,string prefix,string outdir,string sp,int_vector<8>& resultstr){
    
	//string prefix(argv[2]);
	
	vector<string>seqs;
	for(unsigned i = 0;i<dirs.size();i++){
	        string filename = dirs[i];
		ifstream seqfile(filename.c_str());
		string line;
		getline(seqfile,line);
		string seq = "";
		while(seqfile>>line)
		    seq+=line;
		seqfile.close();
		if(seq.length() < 100000)continue;
		seqs.push_back(seq);
		cout<<seq.length()<<endl;
	}
	
        getresultstr(seqs,outdir,sp,resultstr);	
	//gencommon(seqs,outdir,sp);
	
	

	


}
void extract(string parentdir,vector<string>&dirs,unsigned int layer,string taxid){
    if(layer <= 7){
        vector<string>seqs;
        for(unsigned i = 0;i<dirs.size();i++){
            string filename = dirs[i];
            ifstream seqfile(filename.c_str());
            string line;
            getline(seqfile,line);
            string seq = "";
            while(seqfile>>line)
                seq+=line;
            seqfile.close();
            if(seq.length() < 100000)continue;
            seqs.push_back(seq);
            cout<<seq.length()<<endl;
        }
        gencommon(seqs,parentdir,taxid);
    }
 
        
    
}

