#Tyler Kent
#14 March 2017
#Extract degenerate basepair positions

import argparse
import pandas as pd

def arguments():
	parser  = argparse.ArgumentParser(description="Loads chr-pos-sense-bp reference file from bedtools getfasta and returns bed file of degenerate sites. Currently supports only 4fold degeneracy.")
	parser.add_argument("-i", "--input", help="path for input", required=True)
	parser.add_argument("-o","--output", help="path for output", required=True)
	parser.add_argument("-d","--degeneracy", help="Level of degeneracy at which to filter", type=int, default=4)
	args = parser.parse_args()
	return(args)

def reformat(indf):
	
	'''
	split input df into a more usable format
	'''
	
	#split up input to useful dfs
	info = indf[0].apply(lambda x: pd.Series(x.split(":")))
	seq = indf[1]
	chr = info[0]
	start = info[1].apply(lambda x: pd.Series(x.split("-")))[0]
	end = info[1].apply(lambda x: pd.Series(x.split("-")))[1].apply(lambda x: pd.Series(x.split("(")))[0] #lol
	sense = info[1].apply(lambda x: pd.Series(x.split("(")))[1].apply(lambda x: pd.Series(x.rstrip(")")))

	#convert seq to codons
	seq = seq.apply(lambda x: list(x[i:i+3].upper() for i in range(0, len(x), 3)))

	#make lists of 3rd positions conditional on sense and corrected for hanging nucleotides
	possense = pd.concat([sense,start,end],axis=1)
	possense.columns = ['sense','start','end']
	possense['end'] = possense.apply(lambda x: x['end'] if (int(x['end'])-int(x['start'])) % 3 == 0 else int(x['end'])+2 if (int(x['end'])-int(x['start'])) % 3 == 1 else int(x['end'])+1 if (int(x['end'])-int(x['start'])) % 3 == 2 else 'poo',axis=1)
	posseries = possense.apply(lambda x: list(range(int(x['start']),int(x['end']),3))[::-1] if x['sense']=='-' else list(range(int(x['start'])+2,int(x['end']),3)),axis=1)

	#get lists of chromosome name for each position
	chrseq = pd.concat([chr,seq],axis=1)
	chrseq.columns = ['chr','seq']
	chrseries = chrseq.apply(lambda x: [x['chr']]*len(x['seq']),axis=1)

	#stack series to form columns of new dataframe
	seqstacked = seq.apply(lambda x: pd.Series(x)).stack().reset_index(drop=True)
	posstacked = posseries.apply(lambda x: pd.Series(x)).stack().reset_index(drop=True)
	chrstacked = chrseries.apply(lambda x: pd.Series(x)).stack().reset_index(drop=True)

	#create new data frame
	df = pd.concat([chrstacked,posstacked,seqstacked],axis=1)
	df.columns = ['chr','start','codon']
	df['end'] = df['start']+1
	df = df[['chr','start','end','codon']]
	
	return(df)

def is_fourfold(df,degenerate):
	
	'''
	Filter dataframe by fourfold codons
	'''

	filtered = df[df['codon'].isin(degenerate.get('4fold'))]
	bed = filtered[['chr','start','end']]
	bed['start'] = bed['start'].astype(int)
	bed['end'] = bed['end'].astype(int)

	return(bed)

def define_sites():

	'''
	Create look-up dictionary of degenerate codons
	'''

	Degenerate = {
		"4fold":{"TCT","TCC","TCG","TCA","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG",
			"CGT","CGC","CGA","CGG","ACT","ACC","ACA","ACG","GTT","GTC","GTA","GTG",
			"GCT","GCC","GCA","GCG","GGT","GGC","GGA","GGG"
		}
	}
	
	return(Degenerate)


##########################
# Begin
##########################

args = arguments()

#get file info
inpath = args.input
outpath = args.output

Degenerate = define_sites()

chunksize = 500
reader = pd.read_table(inpath, sep='\t', chunksize=chunksize, iterator = True, comment='#', header=None)	

iter=0
for chunk in reader:
	df = reformat(chunk)
	bed = is_fourfold(df,Degenerate)
	if iter==0:
		bed.to_csv(outpath,header=False,index=False,sep='\t')
	else:
		bed.to_csv(outpath,header=False,index=False,mode='a',sep='\t')
	iter+=1
