import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
import os.path
from multiprocessing import Pool

def parsegbk(filename):
    gbfile = open(filename,'r')
  
    for record in SeqIO.parse(gbfile,'genbank'):
	    for f in record.features:
		    if f.type == 'CDS':
			    if 'gene' in f.qualifiers:
				    protid = f.qualifiers['protein_id'][0]
			            geneid = f.qualifiers['gene'][0]
				    print protid+'\t'+geneid+'\t'+str(len(f.location.extract(record)))

    gbfile.close()

def processdir(args,dirname,filenames):
	for filename in filenames:
                info = filename.rstrip().split('.')
                if info[len(info)-1] == 'gbk':
		    print dirname+'/'+filename
		    parsegbk(dirname+'/'+filename)


def job(filename):
	return filename


gbklist = []
for fpath,dirname,filename in os.path.walk(sys.argv[1]):
	info = filename.rstrip().split('.')
	if info[len(info)-1] == 'gbk':
		gbklist.append(dirname+'/'+filename)
print len(gbklist)

p = Pool(6)

results = p.map(job,gbklist)

print len(results)
print results[0]
