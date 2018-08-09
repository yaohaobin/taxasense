import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

gbfile = open(sys.argv[1],'r')
for record in SeqIO.parse(gbfile,'genbank'):
	for f in record.features:
		if f.type == 'CDS':
			if 'gene' in f.qualifiers:
				for protid in f.qualifiers['protein_id']:
					print protid
				print len(f.location.extract(record))

gbfile.close()