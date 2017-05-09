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
	chr = info[2]
	name = info[0]
	start = info[3].apply(lambda x: pd.Series(x.split("-")))[0]
	end = info[3].apply(lambda x: pd.Series(x.split("-")))[1].apply(lambda x: pd.Series(x.split("(")))[0] #lol
	sense = info[3].apply(lambda x: pd.Series(x.split("(")))[1].apply(lambda x: pd.Series(x.rstrip(")")))

	#annoyingly drop duplicate isoforms (keep longest)
#	tmpdf = pd.concat([chr,start,end,name,sense,seq],axis=1,ignore_index=True)
#	tmpdf.columns = ['chr','start','end','name','sense','seq']
#	tmpdf1 = tmpdf.drop_duplicates(subset = 'start',keep='last')
#	tmpdf2 = tmpdf1.drop_duplicates(subset = 'end',keep='last')
#	tmpdf2.reset_index(inplace=True)
#	seq = pd.Series(tmpdf['seq'])
#	sense = pd.Series(tmpdf['sense'])
#	name = pd.Series(tmpdf['name'])
#	end = pd.Series(tmpdf['end'])
#	start = pd.Series(tmpdf['start'])
#	chr = pd.Series(tmpdf['chr'])

	#convert seq to codons
	seq = pd.Series(seq.apply(lambda x: list(x[i:i+3].upper() for i in range(0, len(x), 3))))
	seq.to_csv('checkthis')
	#make lists of 3rd positions conditional on sense and corrected for hanging nucleotides
	possense = pd.concat([sense,start,end],axis=1)
	possense.columns=(['sense','start','end'])
	possense['end'] = possense.apply(endpos,axis=1)
	possense['start'] = possense.apply(startpos,axis=1)
	#possense['end'] = possense.apply(lambda x: x['end'] if (int(x['end'])-int(x['start'])) % 3 == 0 else int(x['end'])+2 if (int(x['end'])-int(x['start'])) % 3 == 1 else int(x['end'])+1 if (int(x['end'])-int(x['start'])) % 3 == 2 else 'poo',axis=1)
#	(possense.apply(lambda x: list(range(int(x['start'])+2,int(x['end']),3)) if x['sense']=='+' else 'poo',axis=1)).to_csv('test')
#	posseries = possense.apply(lambda x: list(range(int(x['start']),int(x['end']),3))[::-1] if x['sense']=='-' else list(range(int(x['start'])+2,int(x['end']),3)) if x['sense']=='+',axis=1)
	try:
#	print(possense)
		posseries = possense.apply(positions,axis=1)
	except ValueError:
		posseries=[]
		for i,item in possense.iterrows():
			posadd = positions(item)
			posseries.append(posadd)	

		posseries = pd.Series(posseries)

	#get lists of chromosome and name for each position
	chrseq = pd.concat([chr,seq,name],axis=1)
	chrseq.columns=(['chr','seq','name'])
	try:
		chrseries = pd.Series(chrseq.apply(lambda x: list([x['chr']]*len(x['seq'])),axis=1))
		nameseries = pd.Series(chrseq.apply(lambda x: list([x['name']]*len(x['seq'])),axis=1))
	except ValueError:
		chrseries=[]
		nameseries=[]
		for i,item in chrseq.iterrows():
			chradd=list([item['chr']]*len(item['seq']))
			nameadd=list([item['name']]*len(item['seq']))

			chrseries.append(chradd)
			nameseries.append(nameadd)
		chrseries = pd.Series(chrseries)
		nameseries = pd.Series(nameseries)

	#stack series to form columns of new dataframe
	seqstacked = seq.apply(lambda x: pd.Series(x)).stack().reset_index(drop=True)
	posstacked = posseries.apply(lambda x: pd.Series(x)).stack().reset_index(drop=True)
	chrstacked = chrseries.apply(lambda x: pd.Series(x)).stack().reset_index(drop=True)
	namestacked = nameseries.apply(lambda x: pd.Series(x)).stack().reset_index(drop=True)

	#create new data frame
	df = pd.concat([chrstacked,posstacked,seqstacked,namestacked],axis=1)
	df.columns = ['chr','start','codon','name']
	df['end'] = df['start']+1
	df = df[['chr','start','end','codon','name']]
	
	return(df)

def drop_alternate_isoforms(df):

	'''
	Convert pandas df to list to iterate over rows
	Check for overlapping rows
	Keep only longest row
	'''

	rows = df.values.tolist()

	iter = 0
#	for row in rows:

		

def endpos(x):

	'''
	Redefine end position if pos sense not divisible by 3
	'''

	if (int(x['end'])-int(x['start'])) % 3 == 0 and x['sense'] == "+":
		return(int(x['end']))

	elif (int(x['end'])-int(x['start'])) % 3 == 1 and x['sense'] == "+":
		return(int(x['end'])+2)

	elif (int(x['end'])-int(x['start'])) % 3 == 2 and x['sense'] == "+":
		return(int(x['end'])+1)

	elif x['sense'] == "-":
		return(int(x['end']))


def startpos(x):

	'''
	Redefine start position if neg sense not divisible by 3
	'''

	if (int(x['end'])-int(x['start'])) % 3 == 0 and x['sense'] == "-":
		return(int(x['start']))

	elif (int(x['end'])-int(x['start'])) % 3 == 1 and x['sense'] == "-":
		return(int(x['start'])-2)

	elif (int(x['end'])-int(x['start'])) % 3 == 2 and x['sense'] == "-":
		return(int(x['start'])-1)

	elif x['sense'] == "+":
		return(int(x['start']))


def positions(x):

	'''
	Get lists of positions from seq data dependent on sense
	'''

	newstart = x['start'] + 2
	newend = x['end'] + 1

	if x['sense'] == "-":
		test = list(range(newstart,newend,3))[::-1]
		test = [x - 2 for x in test]

	elif x['sense'] == "+":
		test = list(range(newstart,newend,3))
#		test = [x + 2 for x in test]
	pd.Series(test).to_csv('checkpos')
	return(test)


def is_fourfold(df,degenerate):
	
	'''
	Filter dataframe by fourfold codons
	'''

	filtered = df[df['codon'].isin(degenerate.get('4fold'))]
	bed = filtered[['chr','start','end','name']]
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

chunksize = 1000
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
