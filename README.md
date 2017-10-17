# Deduper

## Part 1

### Duplicates, UMIs, and soft-clipping

PCR duplicates surface when a transcript becomes amplified during PCR in a weighted manner, not in line with the transcripts actual expression level. If PCR duplicates are not culled the counts of transcripts for which there are PCR duplicates will be falsely represented. If the duplicates are highly "expressed" enough, the differntial expression analysis will be in vain. 


Search for UMIs
    allow some mismatch
    allow for soft mapping?
Count UMIs
    if duplicates, choose the one with the highest quality score

See UMIClusterer in network.py umi-tools
Tom Smith--to do
adjacency dict?


The start postion of a read is considered to be the start of its alignment minus any soft clipped bases. A read aligned at position 500 with cigar 2S98M will be assumed to start at postion 498.

umi_methods.py
```
def find_splice(cigar):
'''Takes a cigar string and finds the first splice position as  an offset from the start. To find the 5' end (read coords) of the junction for a reverse read, pass in the reversed cigar tuple'''

    offset = 0
    # a soft clip at the end of the read is taken as splicing
    # where as a soft clip at the start is not.
    if cigar[0][0] == 4:
        offset = cigar[0][1]
        cigar = cigar[1:]

    for op, bases in cigar:
        if op in (3, 4):
            # N or S: found the splice
            return offset
        elif op in (0, 2, 7, 8):
            # M, D, = or X: reference consuming
            offset += bases
        elif op in (1, 5, 6):
            # I, H, P: non-reference consuming
            continue
        else:
            raise ValueError("Bad Cigar operation: %i" % op)

    return False
```
group.py
```
def get_read_position(read, soft_clip_threshold):
    ''' get the read position (taking account of clipping) '''
    is_spliced = False

    if read.is_reverse:
        pos = read.aend
        if read.cigar[-1][0] == 4:
            pos = pos + read.cigar[-1][1]
        start = read.pos

        if ('N' in read.cigarstring or
            (read.cigar[0][0] == 4 and
             read.cigar[0][1] > soft_clip_threshold)):

            cigar = read.cigar[::-1]
            is_spliced = find_splice(cigar)
    else:
        pos = read.pos
        if read.cigar[0][0] == 4:
            pos = pos - read.cigar[0][1]
        start = pos

        if ('N' in read.cigarstring or
            (read.cigar[-1][0] == 4 and
             read.cigar[-1][1] > soft_clip_threshold)):
            is_spliced = find_splice(read.cigar)

    return start, pos, is_spliced
```
