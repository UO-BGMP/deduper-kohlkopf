#!/usr/local/bin/python3

import argparse
import re
import random
from pprint import pprint

def get_args():
    parser = argparse.ArgumentParser(prog="kinning_deduper.py", description="Removes PCR duplicates from a SAM file")

    parser.add_argument("-f", help="Path to sorted SAM file. <str>", required=True, type=str)
    parser.add_argument("-p", help="If passed as 'True', data is considered paired-ended. <boolean> (def=False)", required=False, type=bool, default=False)
    parser.add_argument("-umi", help="Path to file containing UMI sequences. <str>", required=True, type=str, default='')
    #parser.add_argument("-h", help="Add some helpful information.", required=False, type=str)
    
    #use this in jupyter nb
    #return parser.parse_args("-f Dataset1_sorted.sam -umi umi_list.txt".split())

    #use this in actual script
    return parser.parse_args()
args = get_args()

class Read_Info:
    '''Object with header/read information and a function for adjusting read position to correct for soft clipping.'''
    
    def __init__(self, line):
        self.tag = line.split('\t')[0].split(':')[-1]
        self.flag = int(line.split('\t')[1])
        self.chr = line.split('\t')[2]
        self.pos = int(line.split('\t')[3])
        self.cig = line.split('\t')[5]
        self.seq = line.split('\t')[9]
        self.qual = line.split("\t")[10]
        self.full_line = line

    def adj_pos(self, adj):
        self.pos = self.pos-adj 
        

def umi_list(umi_file):
    '''Returns list of UMIs passed in by -umi'''

    with open(umi_file, 'r') as umis:
        return [line.strip('\n') for line in umis.readlines()]
    
def check_strand(flag, paired):
    '''Checks the bitwise flag and returns the strand of the read. '''

    strand = '+'
    if (flag & 16) == 16: 
        strand = '-'  #if flag indicates read is reverse complemented return "-"

    return strand

def conv_phred(qual):
    '''Converts individual quality scores on a line to phred score, returns the sum.'''

    sum_phreds=0
    
    for letter in qual:
        n = ord(letter) - 33
        sum_phreds+=n
        return sum_phreds

def dup_remover(reads, dups):
    '''Takes list of reads and list of duplicates. Returns a list with duplicates removed. Keeps the version with the highest average Phred score in the batch.'''

    to_keep = []
    best_q = 0
    
    for read in dups:
        read_q = conv_phred(Read_Info(read).qual)
        if read_q > best_q:
            best_q = read_q
            to_keep.append(read)

    dups = [y for y in dups if y not in to_keep]
    kept_reads = [x for x in reads if x not in dups]
    return kept_reads


def dedup_batch(batch, paired):
    '''Takes in batch of reads (list), and uses the first read in the batch as the reference. Uses this read to determine any PCR duplicates of the read in the batch. Each read's POS is adjusted if needed based on the CIGAR field (for soft clipping), and strand orientation is determined from the FLAG before reads are compared. Duplicates are then removed from list of reads, and deduped batch is returned.'''

    #pprint(batch)
    
    kept_reads = batch
    first_read = Read_Info(batch[0])
    dups = []

    for i in range(1, len(batch)):
        pos_dup = Read_Info(batch[i]) #possible duplicate
        if check_strand(first_read.flag, paired) == check_strand(pos_dup.flag, paired): #confirm same strand
            if 'S' in pos_dup.cig and pos_dup.cig.find('S')+1 < len(pos_dup.cig): #check if there is soft-clipping at the 5'-end such that the pos has to be adjusted
                adj = int(pos_dup.cig[:pos_dup.cig.find('S')]) #get value to adjust by from CIGAR
                pos_dup.adj_pos(adj)

            if pos_dup.pos == first_read.pos and pos_dup.tag == first_read.tag: #put in to dups if umi and positions match
                dups.append(pos_dup.full_line)
    
    if len(dups) > 0:
        dups.append(first_read.full_line)
        kept_reads = dup_remover(batch, dups)

    return kept_reads



def read_batch_maker(line, umis, paired, out_file):
    '''Ensures that UMIs are valid. Then lines are read from the file in batches and a single version of the line is stored in list.'''

    first_read = Read_Info(line) # create Read_Info object from FILE line
    
    if umis != []: #if umi file provided
        utag = first_read.tag
        while utag not in umis: #validate umi
            first_read = Read_Info(FILE.readline())
            utag = first_read.tag
        else:
            pass

    r_len = len(first_read.seq) #determine read length and add to first_read pos to find cut-off for position batch
    
    batch = []
    batch.append(first_read.full_line) #add first read

    while first_read: #until eof
        next_read = FILE.readline() #read in next line
        if next_read == '':
            break
        
        next_read = Read_Info(next_read) #create Read_Info object from next_read
        while next_read.pos <= first_read.pos+r_len and next_read.chr == first_read.chr: #check to see if the position of the next read is within the first read's position + read length
            if umis != []: #if umi file provided
                if next_read.tag in umis: #validate umi
                    batch.append(next_read.full_line)
                    next_read = FILE.readline()
                    if next_read == '': #EOF
                        break

                    next_read = Read_Info(next_read)

                else: #if umi invalid, skip this read
                    next_read = FILE.readline()
                    if next_read == '':
                        break

                    next_read = Read_Info(next_read)

        while len(batch) > 1 and first_read.pos+r_len < next_read.pos:
            batch = dedup_batch(batch, paired) #dedup the batch
            with open(out_file,'a') as out:
                out.write(batch[0])
            
            if len(batch) == 1:
                break

            first_read = Read_Info(batch[1])
            batch = batch[1:]
            
######            
            
#         if len(batch) > 1:
#              batch.append(next_read.full_line)
                
#####                
            
        if len(batch) == 1: #once batch is deduplicated, reset the batch
            with open(out_file,'a') as out:
                out.write(batch[0])

            batch = []
            batch.append(next_read.full_line) #begin with the next read, the one just excluded by the batch maker
            first_read = Read_Info(batch[0])




    with open(out_file,'a') as out: 
        out.write(first_read.full_line)

umis = []
if args.umi != '': #if umi file provided, set umis to provided list
    umis = umi_list(args.umi)

with open(args.f, 'r') as FILE:

    filename = re.sub("\/.+\/", "", args.f)
    line = FILE.readline() #read in first line of the SAM file

    header = [] #keep header lines (@)
    while line.startswith('@') == True:
        header.append(line)
        line = FILE.readline()

    if header != []:
        with open(str(filename.strip(".sam"))+'_deduped.sam','a') as out:
            out.write(''.join(header))
            
    output_file = str(filename.strip(".sam"))+'_deduped.sam'
    read_batch_maker(line, umis, args.p, output_file)
