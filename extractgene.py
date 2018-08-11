import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
import os
from multiprocessing import Pool

def parsegbk(filename):
    print filename
    info = filename.split('/')
    lastname = info[len(info)-1]
    #outfile = open(sys.argv[2]+'/'+lastname,'w')

    gbfile = open(filename,'r')
    outstr = filename+'\n'  
    for record in SeqIO.parse(gbfile,'genbank'):
	    for f in record.features:
		    if f.type == 'CDS':
			    if 'gene' in f.qualifiers:
				    protid = f.qualifiers['protein_id'][0]
			            geneid = f.qualifiers['gene'][0]
				    outstr = outstr + protid+'\t'+geneid+'\t'
                                    outstr = outstr + f.location.extract(record).seq+'\n'

    gbfile.close()
    print lastname,len(outstr)
    outfile = open(sys.argv[2]+'/'+lastname,'w')
    outfile.write(str(outstr))
    outfile.flush()
    outfile.close()
    return lastname

def processdir(args,dirname,filenames):
	for filename in filenames:
                info = filename.rstrip().split('.')
                if info[len(info)-1] == 'gbk':
		    print dirname+'/'+filename
		    parsegbk(dirname+'/'+filename)


def job(filename):
           
        return parsegbk(filename)
	
def outputresult(result):
    print result[0]
    outfile = open(sys.argv[2]+'/'+result[0],'w')
    outfile.write(str(result[1]))
    outfile.close()

gbklist = []
for fpath,dirname,filenames in os.walk(sys.argv[1]):
        for filename in filenames:
	    info = filename.rstrip().split('.')
	    if info[len(info)-1] == 'gbk':
		    gbklist.append(fpath+'/'+filename)
print len(gbklist)

p = Pool(12)

idx = 0
p.map_async(job,gbklist)
p.close()
p.join()
#print len(results)


    
    #outfile = open(sys.argv[2]+'/'+lastname,'w')
    #outfile.write(results[idx].get())
    #outfile.close()
   

