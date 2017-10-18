# Deduper

## Part 1

### Duplicates, UMIs, and soft-clipping

PCR duplicates surface when a transcript becomes amplified during PCR in a weighted manner, not in line with the transcripts actual expression level. If PCR duplicates are not culled the counts of transcripts for which there are PCR duplicates will be falsely represented. If the duplicates are highly "expressed" enough, the differential expression analysis will be in vain. 


### Quasicode

#### Functions

```
parse_concise_idiosyncratic_gapped_alignment_report(align_pos, cigar_string):
	if (soft clipping present == TRUE)
		shift align_pos
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
	parse_concise_idiosyncratic_gapped_alignment_report(align_pos, cigar_string)
	if(confirm_umi(umi, umi_dict) == TRUE):
	    check_pair(umi, final_align_pos)
```

