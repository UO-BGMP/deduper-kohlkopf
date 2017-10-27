# Deduper

## Part 1 

### Duplicates, UMIs, and soft-clipping

PCR duplicates in read data surface when a transcript from the same molecule becomes sequenced several times, incorrectly reflecting the actual expression level of the transcript. If PCR duplicates are not culled the counts of transcripts will be false. If the duplicates are highly "expressed" enough, the differential expression analysis will be in vain as we will be comparing artificially expressed transcripts. We can incorporate Unique Molecular Identifiers (UMIs) in to our experimental procedure to identify and cull PCR duplicates in the read data. If we don't have UMIs we can compare the chromosome location, strandedness and alignment position of the reads. Soft clipping occurs when bases are masked in the alignment but remain in the sequence of the SAM file. This shifts the alignment position, causing issues with the the chromosome location, strandedness and alignment position comparison technique. The alignment position simply needs to be shifted to account for the soft clipping


### Quasicode

#### functions()

```
parse_concise_idiosyncratic_gapped_alignment_report(align_pos, cigar_string):
	if (soft clipping present == TRUE)
		final_align_pos = shift align_pos
	else
	    break
	return final_align_pos
	
#ensure actual UMI against list of knowns rather than seq error
confirm_umi(umi, umis):
	if(umi is in umis)
	    status = TRUE
    else
        status = FALSE
	return status
	

check_pair(umi, final_align_pos):
	if(umi:final_align_pos is in dict)
		write formated line to SAM
```

#### main()
```
Create a dictionary with the UMIs as keys, the alignment positions including strandedness associated with them will be the values. This will be a list of UMI alignment pairs.

umi_dict = UMIs and tuple(strand flag, align_pos)

umis = 96 actual UMIs

for line in file:
	final_align_pos = parse_concise_idiosyncratic_gapped_alignment_report(align_pos, cigar_string)
	if(final_align_pos != align_pos):
	    align_pos = final_align_pos
	
	if(confirm_umi(umi, umi_dict) == TRUE):
	    check_pair(umi, final_align_pos)
```

