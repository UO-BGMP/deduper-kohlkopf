import sys
maximumscore=61
import time

#command line parser- homemade to cut down on packages needed in a distro
menu="""

Written by Matthew Burriesci in 2011 as part of a PhD in Genetics.

This program takes fastq data from single or paired end illumina data or 454 reads. It outputs data that has been binned into either files on a local computer or mapped to reducers if you are running mapreduce.

-h print this screen and quit  

-s single machine (default)
-r map reduce

-i [input files with path]
-o [output path] include platform-specific slashes, this is where the output files will be stored

-t [char] read type, 's' for single ended read illumina, 'p' for paired illumina read(default), '4' for 454 data

-m [int] in MB. defaults to 20MB. python get the max size of a single output file. A 140MB file takes 1.9GB of RAM on a 64 bit workstation. The smaller the filesize, the faster a job will run but the less compression will be achieved. Therefore, it is recommended that you take each computer you will be using and take the total available RAM in MB(after operating system footprint has been subtracted) and divide by the number of cpus. then take this number and divide by 20 to get your max file size in MB. You may find that it is sometimes better to remove smaller machines from the cluster or not use the max number of cpus on those macxhines to get the largest file size.

-f [int] first bin for the names of files. 6 is default. not used for map-reduce
-b [int] second bin- the first kmer of this length in each sequence will be used to bin the sequences. only sequences in the same bin can be compared. This is used in map-reduce
-c [int] size of cache in MB. default is 100MB
-u [int] size of read in Bytes. default is 335 Bytes

"""

#main
#parse the command line
seqdict = {}
scoredict = {}
compseq = {}
argsdict ={'-s':0, '-l':0, '-r':0, '-t':'p', '-m':20, '-i':'', '-o':'', '-f':6, '-b':10, '-c':100, '-u':335}
argv=sys.argv
j=0
while (j<len(argv)-1):
	j+=1
	if argv[j] == '-h':
		print menu
		sys.exit(-1)
	elif argv[j] == '-s':
		if argsdict['-r']!=0:
			print "too many single/mapreduce args"
			sys.exit(-1)
		argsdict['-s']=1
	elif argv[j] == '-l':
		
		print "no local network setting for local parser"
		sys.exit(-1)
	elif argv[j] == '-r':
		if argsdict['-s']!=0 or argsdict['-l']!=0:
			print "too many single/local network/mapreduce args"
			sys.exit(-1) 
		else:
			argsdict['-r']=1
	elif argv[j] == '-t':
		if argv[j+1]=='s' or argv[j+1]=='p' or argv[j+1]=='4':
			argsdict['-t']=argv[j+1]
			j+=1
		else:
			print "invalid read type specified"
			sys.exit(-1)
	elif argv[j] == '-i':
		argsdict['-i']=argv[j+1]
		j+=1
	elif argv[j] == '-o':
		argsdict['-o']=argv[j+1]
		j+=1
	elif argv[j] == '-m':
		argsdict['-m']=argv[j+1]
		j+=1
	elif argv[j] == '-f':
		argsdict['-f']=argv[j+1]
		j+=1
	elif argv[j] == '-b':
		argsdict['-b']=argv[j+1]
		j+=1
	elif argv[j] == '-c':
		argsdict['-c']=argv[j+1]
		j+=1
	elif argv[j] == '-u':
		argsdict['-u']=argv[j+1]
		j+=1
	else:
		print "unknown identifier in commandline:"
		print argv[j]
		sys.exit(-1)
	
if (argsdict['-o']=='' or argsdict['-i']=='') and argsdict['-r']!=1:
	print "please define in file and out folder"
	sys.exit(-1)

#paired end illumina functions
def compressor(templine, kmer_len, kmer_len2):
    num=0
    bins[templine[0][0:int(kmer_len)]].append(templine[0][0:int(kmer_len2)]+"\t%0.0f:%s:%s:%s:%s" %(x,templine[0], templine[1], templine[2], templine[3]))

#prints data into the appropriate bins, during the first pass (datafile>>many_smaller_bins)
def printer(outfile):
	print time.clock()
	for bin in bins:
		openprint=open(outfile+bin, 'a')
		if bins[bin] !=list():
			for seq in bins[bin]:
				print >> openprint, seq
			openprint.close()
			bins[bin]=list()
#prints the contents of a bin file into smaller buckets when necessary due to large size of original bin
def printer2(outfile, bigbin):
	print time.clock()
	for bin in bins:
		if bins[bin] !=list():
			if len(bins[bin])>100:
				openprint=open(outfile+bin, 'a')
				for seq in bins[bin]:
					print >> openprint, bin, "\t", seq
				openprint.close()
				bins[bin]=list()
			else:
				openprint=open(outfile+bigbin, 'a')
				for seq in bins[bin]:
					print >>openprint, bin, "\t", seq
				openprint.close()
				bins[bin]=list()

#single end illumina/454 functions
def compressorsingle(templine, kmer_len, kmer_len2):
    num=0
    bins[templine[0][0:int(kmer_len)]].append(templine[0][0:int(kmer_len2)]+"\t%0.0f:%s:%s" %(x,templine[0], templine[1]))

#sets some variables from the command line args
size_of_cache=int(argsdict['-c'])#per processor, in MB
kmer_len=int(argsdict['-f'])
kmer_len2=int(argsdict['-b'])

#build cache of processed reads. create bins, each bin being a kmer of length specified in arg. 10 bases gives 9.7M bins (~1M non N containing)
bins = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
for m in range(1,int(kmer_len)):
	tempbin={}
	for bin in bins:
		current_key=bin
		for base in ['A', 'T', 'G', 'C', 'N']:
			tempbin[current_key+base] =list()
	bins = tempbin
		



#needed only for non map-reduce work
	import dircache
	import os
	file1 =open(argsdict['-i'] ,"r")
	indexingfile = open("index.txt","w")

	if argsdict['-t']=='p': #paired end illumina
		#process file, merging 8? lines into 1 tab and colon seperated file
		x =0
		z=0
		temp_line=[None]*4
		indexing=''
		space=0
		for line in file1:
		    line = line.strip()
		    if (line[0] == '@' and line.find('#0/1')!=-1):#print a completed line
			#print line
			if x>0:
			    compressor(temp_line, kmer_len, kmer_len2)
			    print >>indexingfile, "%0.0f\t%s" %(x,indexing)
			x=x+1
			z+=1
			temp_line=[None]*4
			space=0
			splitline=line.split('#')
			indexing=splitline[0]
		    if (line[0] != '@' and line[0] !='+'):
			temp_line[space] = line
			space=space+1
			if (z*int(argsdict['-u']))>1000000*int(size_of_cache):
				z=0
				print str(x) +" read processed"#prints status for user
				printer(argsdict['-o'])	
		#print out last line
		compressor(temp_line, kmer_len, kmer_len2)
		printer(argsdict['-o'])
		print >>indexingfile, "%0.0f\t%s" %(x,indexing)
		file1.close()
	else:#single end or 454 first processing	
		#process file 4 lines of fastq into 1 line
		x =0
		z=0
		temp_line=[None]*2
		indexing=''
		space=0
		firstline=1
		for line in file1:
		    line = line.strip()
		    if (line[0] == '@'):#print a completed set of lines
			#print line
			if x>0:
			    compressorsingle(temp_line, kmer_len, kmer_len2)
			    print >>indexingfile, "%0.0f\t%s" %(x,indexing)
			x=x+1
			z+=1
			temp_line=[None]*2
			space=0
			indexing=line
		    if (line[0] != '@' and line[0] !='+'):
			#if space=0: #trim 
			temp_line[space] = line
			space=space+1
			if (z*int(argsdict['-u']))>1000000*int(size_of_cache):
				z=0
				print str(x) +" read processed"
				printer(argsdict['-o'])	
		#print out last lines
		compressorsingle(temp_line, kmer_len, kmer_len2)
		printer(argsdict['-o'])
		print >>indexingfile, "%0.0f\t%s" %(x,indexing)
		file1.close()
	
	#start processing files to make them smaller, splitting up larger files into more manageable sizes by taking sequences with different kmers in the same file and splitting out the most overrepresented
	inputs=dircache.listdir(argsdict['-o'])#get all files
	for input in inputs:
		size=os.path.getsize(argsdict['-o']+input)
		if int(argsdict['-m'])*1000000<size:#file too big
			file2 =open(argsdict['-o']+input,"r")
			x=0
			z=0
			bins={}
			for line in file2:
				x+=1
				z+=1
				line = line.strip()
				key=line.split('\t')#key is the line (why split on \t?)
				if 0==bins.has_key(key[0]):#if key does not exist create bin and add it
					bins[key[0]] = []#initialize new bin within bins
					bins[key[0]].append(key[1])#add key[1] to the the bin
	    			else:
					bins[key[0]].append(key[1])#if key does exist add key[1]
				if (z*int(argsdict['-u']))>1000000*int(size_of_cache):
					z=0
					print str(x) +" reads put into smaller files "+input
					printer2(argsdict['-o'], input+'a')
			print str(x) +" reads put into smaller files "+input
			printer2(argsdict['-o'], input+'a')	
			file2.close()
			os.remove(argsdict['-o']+input)
		if size==0:#empty, so delete it
			os.remove(argsdict['-o']+input)
	#process files with same large bins into smaller files (buckets)- this should only happen with a couple of bins and can reduce the final amount of compression
	inputs=dircache.listdir(argsdict['-o'])#get all files again
	for input in inputs:
		size=os.path.getsize(argsdict['-o']+input)
		if int(argsdict['-m'])*1000000<size:#file too big
			file2 =open(argsdict['-o']+input,"r")
			x=0
			z=0
			save=[]
			suffix=0
			for line in file2:
				x+=1
				z+=1
				line= line.strip()
				save.append(line)
				if ((z*int(argsdict['-u']))>1000000*int(size_of_cache) or (z*int(argsdict['-u']))>1000000*int(argsdict['-m'])):
					z=0
					outputfile=open(argsdict['-o']+input+str(suffix), 'w')
					for saveline in save:
						print >>outputfile, saveline
					print str(x) +" reads put into smaller files "+input
					outputfile.close()
					suffix+=1
					save=[]
			outputfile=open(argsdict['-o']+input+str(suffix), 'w')
			for saveline in save:
				print >>outputfile, saveline
			print str(x) +" reads put into smaller files "+input
			outputfile.close()
			suffix+=1
			save=[]
			file2.close()
			os.remove(argsdict['-o']+input)
