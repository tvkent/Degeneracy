from collections import OrderedDict
import argparse
import pandas as pd

def arguments():
	parser  = argparse.ArgumentParser(description="Filters bedtools getfasta tab file to only keep longest isoforms")
	parser.add_argument("-i", "--input", help="path for input", required=True)
	parser.add_argument("-o","--output", help="path for output", required=True)
	args = parser.parse_args()
	return(args)

def order_check(lines,outpath):

	'''
	Loop over lines and check if next line's start position 
	overlaps with previous line's span.
	Send overlapping sections or single nonoverlapping lines
	to function to choose the longest line to write to file.
	'''

#	overlaps = OrderedDict()
#	overlaps = pd.DataFrame(columns = ['line','start','end','span'])
	last = lines[-1]
	fniter = 0
	iter = 0
	
	for line in lines:
		#set intial line to previous line info for the first comparison
		if iter==0:
			sline = line.split()[0].split(":")
			pchr = sline[2]
			pstart = sline[3].split("-")[0]
			pend = sline[3].split("-")[1].split("(")[0]
			pspan = int(pend)-int(pstart)
			overlaps = pd.DataFrame(columns = ['line','start','end','span'])
			overlaps = overlaps.append(pd.Series([line,pstart,pend,pspan],index=['line','start','end','span']),ignore_index=True)
#			check=False
#			newkey = "-".join(pstart,pend,pspan)
#			overlaps[newkey]  = line

		else:
			chr = line.split()[0].split(":")[2]
			
			#only split and compare things if they're on the same chromosome
			if chr == pchr:
				sline = line.split()[0].split(":")
				start = sline[3].split("-")[0]
				end = sline[3].split("-")[1].split("(")[0]
				
				if int(start) < int(pend):
					span = int(end)-int(start)
					lineseries = pd.Series([line,start,end,span],index=['line','start','end','span'])
					overlaps = overlaps.append(lineseries,ignore_index=True)
#					overlaps[newkey] = line
					
					if line == last:
						keep_longest(overlaps,outpath,fniter)
#						fniter+=1
#						check=True

				else:
					keep_longest(overlaps,outpath,fniter)
					fniter+=1
#					check=True

					sline = line.split()[0].split(":")
					pchr = sline[2]
					pstart = sline[3].split("-")[0]
					pend = sline[3].split("-")[1].split("(")[0]
					pspan = int(pend)-int(start)
					overlaps = pd.DataFrame(columns = ['line','start','end','span'])
					overlaps = overlaps.append(pd.Series([line,pstart,pend,pspan],index=['line','start','end','span']),ignore_index=True)
#					lineseries = pd.Series([line,pstart,pend,pspan],index=['line','start','end','span'])
#					overlaps = overlaps.append(lineseries,ignore_index=True)

					if line == last:
						outfile = open(outpath,	'a')
						outfile.write(line)
						outfile.close()

			else:
				keep_longest(overlaps,outpath,fniter)
				fniter+=1
#				check=True

				sline = line.split()[0].split(":")
				pchr = sline[2]
				pstart = sline[3].split("-")[0]
				pend = sline[3].split("-")[1].split("(")[0]
				pspan = int(pend)-int(start)
				overlaps = pd.DataFrame(columns = ['line','start','end','span'])
				overlaps = overlaps.append(pd.Series([line,pstart,pend,pspan],index=['line','start','end','span']),ignore_index=True)
#				lineseries = pd.Series([line,pstart,pend,pspan],index=['line','start','end','span'])
#				overlaps = overlaps.append(lineseries,ignore_index=True)


				if line == last:
						outfile = open(outpath,	'a')
						outfile.write(line)
						outfile.close()

		iter+=1


def keep_longest(overlaps,outpath,iter):

	'''
	Take dictionary of span lengths, get the longest one,
	print the corresponding line to the output file.
	Append to file unless the iteration is the first.
	'''
	overlaps['span'] = overlaps['span'].astype(float)
	df = overlaps.sort_values(['span'], ascending=False)
	line = df['line'].iloc[0]
#	keys = list(overlaps.keys())
#	intkeys = list(map(int, keys))
#	keylist = [x.split("-")[2]

#	longest = max(intkeys)

#	longestline = overlaps[str(longest)]

	if iter == 0:
		outfile = open(outpath, 'w')
		outfile.write(line)

	else:
		outfile = open(outpath, 'a')
		outfile.write(line)

	outfile.close()

#####################
# BEGIN
#####################

args = arguments()

inpath = args.input
outpath = args.output

infile = open(inpath, 'r')
lines = infile.readlines()

#pass all lines to function
order_check(lines,outpath)
