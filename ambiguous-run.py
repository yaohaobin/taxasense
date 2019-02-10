import pandas as pd
import sys
import threading
import scipy.stats as ss
import bisect
#code from https://www.geeksforgeeks.org/binary-search-tree-set-1-search-and-insertion/
class Node: 
    def __init__(self,key): 
        self.left = None
        self.right = None
        self.val = key 
        self.leftsize = 1

def insert(root,node): 
    if root is None: 
        root = node 
    else: 
        if root.val < node.val: 
            if root.right is None: 
                root.right = node 
            else: 
                insert(root.right, node) 
        else: 
            root.leftsize += 1
            if root.left is None: 
                root.left = node
                    
            else: 
                insert(root.left, node) 
    
def getrank(root,x):
    if root is None:
        return 0
    if root.val == x:
        return root.leftsize
    elif root.val < x:
        return root.leftsize + getrank(root.right,x)
    else:
        return getrank(root.left,x)
                


class sort_predata (threading.Thread):
    def __init__(self,aroid):
        threading.Thread.__init__(self)
        self.aroid = aroid
    def run(self):
        #threadLock.acquire()
        
        prehit = {}

        predata = aro_pre[self.aroid]
        mixlist = []
        threshold = aro2param[self.aroid]
        totalnum = 0

        for preid,data in predata.items():
            for point in data:
                totalnum += 1
                if  True or point[0] > 94 or point[1]>= threshold:
                    mixlist.append((preid,point[0],point[1]))
        mixlist = sorted(mixlist,key=lambda x:x[1])
        
      
        truenum = 0        
        root = None
        ambisum = 0
        ambinum = 0
        minscore = 10000
        scorelist = []

        aroinother[self.aroid] = {}
       


        for i in range(len(mixlist)):
            preid = mixlist[i][0]
            if mixlist[i][0] == self.aroid:
                truenum += 1
                '''
                if root is None:
                    root = Node(mixlist[i][2])
                else:
                    insert(root,Node(mixlist[i][2]))
                '''
                bisect.insort_right(scorelist,mixlist[i][2])
                if mixlist[i][2] < minscore:
                    minscore = mixlist[i][2]
            else:
                '''
                if not root is None:
                    num = getrank(root,mixlist[i][2])
                    if num>0:
                        ambisum += num
                        ambinum += 1
                '''
                
                if mixlist[i][2] > minscore:
                    
                    ambinum += 1
                    #print(mixlist[i])                    
                    pos = bisect.bisect_right(scorelist,mixlist[i][2])
                    if 0<pos and pos<=len(scorelist):
                        ambisum += pos


                    #threadLock.acquire()
                    '''
                    if preid in aroinother.keys():
                        if self.aroid in aroinother[preid].keys():
                            aroinother[preid][self.aroid] += pos
                        else:
                            aroinother[preid][self.aroid] = pos
                    else:
                        aroinother[preid] = {}
                        aroinother[preid][self.aroid] = pos
                    
                    if preid in aroinother[self.aroid].keys():
                        aroinother[self.aroid] += pos
                    else:
                        aroinother[self.aroid] = pos
                    '''
                    #threadLock.release()

                        
        threadLock.acquire()

        print(self.aroid,len(mixlist),truenum,ambisum,ambinum,totalnum)
        threadLock.release()



prevalence = pd.read_csv(sys.argv[1],header=None,sep='\t')
aroidx = pd.read_csv(sys.argv[2],header=0,sep='\t')

aroparam = pd.read_csv(sys.argv[3],header=None,sep='\t')
aro2param = {}
for idx,row in aroparam.iterrows():
    aro2param['ARO:'+str(int(row[0])).strip()] = row[1]


aroidlist = []
id2name = {}
aroinother = {}
for idx,row in aroidx.iterrows():
    aroidlist.append(row['ARO Accession']) 
    id2name[row['ARO Accession']] = row['ARO Name']

aro_pre = {}

nullnum = 0
for idx,row in prevalence.iterrows():
    aroid = row[1].split('|')[2]
    
    preinfo = row[0].split('|')
    if len(preinfo) < 3:
        nullnum += 1
        continue
    preid = preinfo[2]
    pid = float(row[2])
     
    score = float(row[11])
    
    point = (pid,score)
    if aroid in aro_pre.keys():
        predata = aro_pre[aroid]
        if preid in predata.keys():
            predata[preid].append(point)
        else:
            predata[preid] = [point]
    else:
        predata = {}
        predata[preid] = [point]
        aro_pre[aroid] = predata


threadLock = threading.Lock()



threadnum = int(sys.argv[4])
i = 0

while i<len(id2name):
     threads  = []
     for j in range(threadnum):
         aroid = aroidlist[i]
         if aroid in aro_pre.keys():
             #print(aroid)
             
             arothread = sort_predata(aroid)
             arothread.start()
             threads.append(arothread)
           
 
         i += 1        
         if i>= len(id2name):
             break

     
     for t in threads:
         t.join()

#print('exit')

prehitout = open(sys.argv[5],'w')
for preid in aroinother.keys():
     outstr = ''
     outstr += preid+'\t'+str(len(aroinother[preid]))
     for aroid,v in aroinother[preid].items():
         outstr += '\t'+aroid+'\t'+str(v)
     prehitout.write(outstr+'\n')

prehitout.close()
