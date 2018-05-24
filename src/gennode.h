#include<vector>
#include<map>
#include<set>
#include<stack>
#include<algorithm>
#include<omp.h>

#include<sdsl/suffix_trees.hpp>
#include"extractcommon.h"
using namespace std;
using namespace sdsl;
/*
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


*/

string findcommon_seq(vector<string> dirs){
	
	cout<<dirs.size()<<endl;
	return ""+(char)dirs.size();
}
string findcommon_node(vector<int>& children,vector<string >&commonstring){
	
	cout<<"common: "<<children.size()<<endl;
	
	return ""+(char)children.size();
}
bool compair(pair<int,int>& p1,pair<int,int>&p2){
    return p1.first<p2.first;
}

void merge_interval(vector<pair<int,int> >&intervals,vector<pair<int,int> >&merged){
   
     if(intervals.size() <= 1)return;
     sort(intervals.begin(),intervals.end(),compair);
     stack<pair<int,int> > stack_merge;
     stack_merge.push(intervals[0]);
     for(unsigned i = 1;i<intervals.size();i++){
         pair<int,int>top_pair = stack_merge.top();
         if(stack_merge.top().second >= intervals[i].first){
             stack_merge.pop();
             stack_merge.push( make_pair(top_pair.first,intervals[i].second) );

         }
         else
            stack_merge.push( intervals[i]);
     }

     merged.resize(stack_merge.size());
     for(unsigned i= 0;i<stack_merge.size();i++){
         merged[i] = stack_merge.top();
         stack_merge.pop();
     }
}
string get_seq(string dir){
	
	return dir; 

}


class commonnode{
public:
   	commonnode* parent;
	
	vector<commonnode*> children;
	
	int id;
	
	bool isLeaf;
	
	int leafnum;
	
	string dir;
        string gi;
	
};

class pathnode{
public:
    vector<int>path;
	

	
	bool isLeaf = false;
	
	int id;
	
};
bool compnode(commonnode* &n1,commonnode* &n2){
        return n1->leafnum > n2->leafnum;
}


class pathtree{
public:
    vector<pathnode*> paths;
	pathnode* root;
	
	
};

class Subphytree{
public:
	vector<commonnode*> commontree; 
	vector<string> commonstring;
	pathtree heavy;
	map<int,string> id_name;

   
	
	commonnode* root;
	int rootid;
	
	void genTree(vector<map<string,set<string> > >& taxtree,map<string,string>& gbkdir);
	
    void heavyPath();
	
	void preorder();
	
	void sortchild();
	void common(string prefix);	
        void commonsp(string prefix,string dbfile,string outdir);
	~Subphytree(){
		for(int i=0;i<commontree.size();i++)
			delete commontree[i];
		for(int i=0;i<heavy.paths.size();i++)
			delete heavy.paths[i];
	}
	
private:
   int label(commonnode* node,int idx);
   void constructHeavypath(commonnode* node,pathnode* newpath,int idx);
};

//construct tree where each node presents commonstrings among subtrees,return the number of nodes (root id = number of nodes - 1)
void Subphytree::genTree(vector<map<string,set<string> > >& taxtree,map<string,string>& gbkdir){

	
	
	
	
	
    
        
    map<string,commonnode*>tempmap;
	
    for(int i=taxtree.size()-1;i>=0;i--){
        for(map<string,set<string> >::iterator itr = taxtree[i].begin();itr!=taxtree[i].end();itr++){
          
			
			commonnode* parentnode = new commonnode;
			parentnode->isLeaf = false;
			parentnode->leafnum = 0;
			tempmap[itr->first] = parentnode;
		        parentnode->dir = itr->first;	
			for(set<string>::iterator itr2=itr->second.begin();itr2!=itr->second.end();itr2++){
				if(gbkdir.find(*itr2) != gbkdir.end() ){
                                         cout<<gbkdir[*itr2]<<endl;
					 commonnode* node = new commonnode;
					 node->isLeaf = true;
					 node->leafnum = 1;
					 node->dir = gbkdir[*itr2];
					 node->gi = *itr2;
					 
					 parentnode->children.push_back(node);
					 node->parent = parentnode;
					 parentnode->leafnum++;
					 
					 
					 
				}
				else{
					 commonnode* node = tempmap[*itr2];
					 parentnode->children.push_back(node);
					 parentnode->leafnum += node->leafnum;
					 node->parent = parentnode;
				}
			} 
			
			
            
        }
    }
    
	root = tempmap[taxtree[0].begin()->first];
        cout<<tempmap.size()<<endl; 
        	
	preorder();
	
	for(map<string,commonnode*>::iterator itr=tempmap.begin();itr!=tempmap.end();itr++){
		id_name[itr->second->id] = itr->first;
	}
}




void Subphytree::sortchild(){
	for(int i=0;i<commontree.size();i++)
		sort(commontree[i]->children.begin(),commontree[i]->children.end(),compnode);
}

void Subphytree::preorder(){
	
        	
	label(root,0);
	sortchild();

        heavyPath();
}

int Subphytree::label(commonnode* node,int idx){
	//cout<<idx<<" "<<node->children.size()<<" "<<node->dir<<endl;
	for(int i=0;i<node->children.size();i++){
		commonnode* child = node->children[i];
                
	        if(child->isLeaf){
                        //cout<<i<<endl;
		        commontree.push_back(child);
			child->id = idx;
			idx++;
		}
		else{
                        //cout<<i<<endl;
			idx = label(child,idx);
			
			idx++;
		}							
	}
	node->id = idx;
	commontree.push_back(node);
	
	return idx;
	
}


void Subphytree::constructHeavypath(commonnode* node,pathnode* newpath,int idx){
	
	
	
	newpath->path.push_back(node->id);
	if(node->isLeaf){
            //for(unsigned i=0;i<newpath->path.size();i++)
                //cout<<newpath->path[i]<<" ";
            //cout<<endl;        
            return;
	}
	bool heavied = false;
	for(int i=0;i<node->children.size();i++){
		
		if(!heavied){
			heavied = true;
			constructHeavypath(node->children[i],newpath,idx);
		}
		else{
			idx++;
			
			pathnode* anotherpath = new pathnode;
			anotherpath->id = idx;
			heavy.paths.push_back(anotherpath);
                        
			constructHeavypath(node->children[i],anotherpath,idx);
		}
	}
}
void Subphytree::heavyPath(){
	
	pathnode* newpath = new pathnode;
	newpath->id = 0;
	heavy.paths.push_back(newpath);
	
	constructHeavypath(root,newpath,0);
	
	
		
}


void Subphytree::common(string prefix){
       vector<int>spnode;
       for(unsigned int i=0;i<commontree.size();i++)
	       for(unsigned int j=0;j<commontree[i]->children.size();j++)
		       if(commontree[i]->children[j]->isLeaf){
			       spnode.push_back(i);
			       break;
		       }
       cout<<"sp num: "<<spnode.size()<<endl;
       for(unsigned int i=0;i<spnode.size();i++){
	       if(commontree[spnode[i]]->children.size() >= 2){
		       cout<<commontree[spnode[i]]->id<<" "<<commontree[spnode[i]]->children.size()<<endl;
                       vector<string>dirs;
                       for(unsigned int j=0;j<commontree[spnode[i]]->children.size();j++)
                           dirs.push_back(commontree[spnode[i]]->children[j]->dir);
                       //extract(dirs,prefix);
               }
       }
       
	
}
void Subphytree::commonsp(string prefix,string dbfile,string outdir){
      ifstream dbtax(dbfile.c_str());
      string gi,ti,sp,ge,fa,ord;
      map<string,string>gisp;
      while(dbtax>>gi>>ti>>sp>>ge>>fa>>ord){
          gisp[gi] = sp;
      }
      dbtax.close();

      map<string,vector<string> > refspdir;
      for(unsigned int i=0;i<commontree.size();i++)
               for(unsigned int j=0;j<commontree[i]->children.size();j++)
                       if(commontree[i]->children[j]->isLeaf){
                               string gi = commontree[i]->children[j]->gi;
                               if(gisp.find(gi) == gisp.end())
                                   continue;
                               string sp = gisp[gi];
                               refspdir[sp].push_back(commontree[i]->children[j]->dir);
                               
                       }
      cout<<"sp num: "<<refspdir.size()<<endl;
      int num = 0;
      
      vector<pair<string,vector<string> > > spdirforomp;
 
      for(map<string,vector<string> >::iterator itr = refspdir.begin();itr!=refspdir.end();itr++){
          spdirforomp.push_back(make_pair(itr->first,itr->second));
          num += itr->second.size();
      }
      cout<<num<<endl;
      
      

      csa_bitcompressed<>csa;
      csa_wt<>csawt;
      //load_from_file(csa,"gesp.csa");
      

      int_vector<16>belong;

      ofstream spseq("spseq.txt");
      
      vector<string>spnums(spdirforomp.size());
      vector<int_vector<8> >results(spdirforomp.size()); 
      #pragma omp parallel for
      for(unsigned int i = 0;i<spdirforomp.size();i++){
          if(spdirforomp[i].second.size()<2)continue;
         
          cout<<"sp "<<spdirforomp[i].first<<" "<<spdirforomp[i].second.size()<<endl;
          vector<string>dirs;
          for(unsigned int j=0;j<spdirforomp[i].second.size();j++)
              dirs.push_back(spdirforomp[i].second[j]);
          
          
          extract(dirs,prefix,outdir,spdirforomp[i].first,results[i]);
          if(results[i].size() !=0)spnums[i] = spdirforomp[i].first;
      }
      cout<<"sum up"<<endl;      
      //ofstream gespcommon("gesp.index");
      int_vector<8> seqtext;
      uint64_t oldlen = 0;
      uint64_t idx = 0;
      uint16_t clade = 0;
      for(unsigned int i=0;i<results.size();i++){
          if(results[i].size() == 0)
              continue;
          spseq<<spnums[i]<<endl;
          oldlen = seqtext.size();
          seqtext.resize(seqtext.size()+results[i].size());
          belong.resize(seqtext.size());
          for(idx=oldlen;idx<seqtext.size();idx++){
              belong[idx] = clade;
              seqtext[idx] = results[i][idx-oldlen];
          }
          clade++;
          
          cout<<seqtext.size()<<endl;
              
          
      }

      
      spseq.close();
      /*(
      vector<uint16_t>table(csa.size());
      table[0] = 0;

      for(uint64_t i = 1;i<csa.size();i++)
          table[i] = belong[csa[i]];  

         
      ofstream fout("tablege");
      for(uint64_t i = 1;i<csa.size();i++)
          fout<<table[i]<<" ";
      fout<<endl;
      fout.close();
      */
      cout<<"wt"<<endl;
      csa_wt<>wt;

      
      construct_im(wt,seqtext);
      store_to_file(wt,"sp.csawt");

      cout<<"csa"<<endl;
      //csa_bitcompressed<>csa;
      construct_im(csa,seqtext);
      store_to_file(csa,"sp.csa");
      
      vector<uint16_t>table(csa.size());
      table[0] = 0;

      for(uint64_t i = 1;i<csa.size();i++)
          table[i] = belong[csa[i]];


      ofstream fout("tablespt");
      for(uint64_t i = 1;i<csa.size();i++)
          fout<<table[i]<<" ";
      fout<<endl;
      fout.close();
      
      
}
