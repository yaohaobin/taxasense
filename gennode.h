#include<vector>
#include<map>
#include<set>
#include<stack>
#include<algorithm>
#include<omp.h>
#include"extractcommon.h"

#define layer 3
using namespace std;

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

class nodeinfo{
public:
	string parentdir;
	vector<string>childrendir;
	string taxid;


};
class commonnode{
public:
   	commonnode* parent;
	
	vector<commonnode*> children;
	
	int id;
	
	bool isLeaf;
	
	int leafnum;
	
	string dir;
	string commonstr;
	
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
	vector<set<string> >taxonomy;
	pathtree heavy;
	map<int,string> id_name;
    map<string,int> tax_id;
   
	
	commonnode* root;
	int rootid;
	
	void genTree(vector<map<string,set<string> > >& taxtree,map<string,string>& gbkdir);
    
	void tree_common();
	void init_common(map<string,string>& gidir,map<string,vector<string> >&dbtax,unsigned int maxlevel,string dbprefix,string commonprefix );
    void heavyPath();
	
	void preorder();
	
	void sortchild();

	void common(string prefix);	

	void readdb(string dbfile,string taxfile,map<string,string>& gidir,map<string,vector<string> >& dbtax,int maxlevel);
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

        //heavyPath();
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
			node->leafnum++;
		}
		else{
                        //cout<<i<<endl;
			idx = label(child,idx);
			node->leafnum+=child->leafnum;
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

void Subphytree::readdb(string dbfile,string taxfile,map<string,string>& gidir,map<string,vector<string> >& dbtax,int maxlevel){
       ifstream db(dbfile.c_str());
       ifstream tax(taxfile.c_str());

    
       


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


void Subphytree::init_common(map<string,string>& gidir,map<string,vector<string> >&dbtax,unsigned int maxlevel,string dbprefix,string commonprefix ){

       //init taxonomy
	   for(unsigned int i=0;i<maxlevel;i++){
	   	   set<string> emptylayer;
	   	   taxonomy.push_back(emptylayer);
	   }


	   //construct leaf (sequence)
	   map<string,commonnode*>tempmap;
       for(map<string,string>::iterator itr=gidir.begin();itr!=gidir.end();itr++){
       	   commonnode* leafnode = new commonnode;
           leafnode->isLeaf = true;
           leafnode->leafnum = 1;
           leafnode->dir=itr->second;
           tempmap[itr->first] = leafnode;

       }

       cout<<"temp size: "<<tempmap.size()<<endl;
       //construct internal nodes (taxonomy)
       for(map<string,vector<string> >::iterator itr = dbtax.begin();itr!=dbtax.end();itr++){
       	   

       	   commonnode* childnode = tempmap[itr->first];
           
       	   for(unsigned int level=0;level<maxlevel;level++){
               
       	   	   if (itr->second[level] == "NA")
       	   	   	   continue;
       	   	   taxonomy[level].insert(itr->second[level]);
       	   	   map<string,commonnode*>::iterator itr2 = tempmap.find(itr->second[level]);
       	   	   commonnode* parentnode;
       	   	   if(itr2==tempmap.end()){
       	   	   	  
       	   	   	   parentnode = new commonnode;
       	   	   	   tempmap[itr->second[level]] = parentnode;
       	   	   	   parentnode->children.push_back(childnode);
       	   	   	   parentnode->isLeaf = false;
       	   	   	   //parentnode->leafnum=1;
       	   	       parentnode->dir =  commonprefix+"/"+itr->second[level]+".common";
       	   	   }
       	   	   else{
       	   	   	   
       	   	   	   parentnode = itr2->second;

                 
       	   	   	   parentnode->children.push_back(childnode);
       	   	   	   
       	   	   	   break;
       	   	   	   //parentnode->leafnum+=childnode->leafnum;
       	   	   }
       	   	   childnode = parentnode;
       	   }
       }
       root = tempmap["1"];
       cout<<"db size: "<<tempmap.size()<<endl;
       preorder();
       for(map<string,commonnode*>::iterator itr=tempmap.begin();itr!=tempmap.end();itr++){
		   id_name[itr->second->id] = itr->first;
		   tax_id[itr->first] = itr->second->id;
	   } 
	   tree_common();
}


void Subphytree::tree_common(){
       for(unsigned int i=0;i<1;i++){
            vector<nodeinfo> layerblock;
       
       	    for(set<string>::iterator itr = taxonomy[i].begin();itr!=taxonomy[i].end();itr++){

       	    	commonnode* taxnode = commontree[tax_id[*itr]];
       	    	nodeinfo info;
                info.parentdir = taxnode->dir;
                info.taxid = *itr;
       	    	for(unsigned int child = 0;child<taxnode->children.size();child++){
       	    		cout<<taxnode->children[child]->dir<<endl;
       	    		info.childrendir.push_back(taxnode->children[child]->dir);
       	    	}
       	    	if(info.childrendir.size() >= 2){
       	    	    
       	    	    layerblock.push_back(info);
       	    	}
       	    	//extract(taxnode->dir,dirs,prefix);
       	    }
       	    #pragma omp parallel for
       	    for(unsigned int j=0;j<layerblock.size();j++){
       	    	cout<<layerblock[j].parentdir<<" "<<layerblock[j].childrendir.size()<<endl;

       	    	extract(layerblock[j].parentdir,layerblock[j].childrendir,i,layerblock[j].taxid);
       	    }
       }
}




