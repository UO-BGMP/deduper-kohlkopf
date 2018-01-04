def strand_checker(bit):
	'''Takes the bitwise flag and checks it for strandedness. Assumes read is mapped and data are single-stranded. Returns "+" or "-", depending on the strand.'''

	if (bit & 4) != 4:
		#raise NameError("Read is unmapped!")	
		return None

	strand = "+"
	if (bit & 16) == True:
		strand = "-"

	return strand

print(strand_checker(4))
