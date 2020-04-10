import sys

#usage: primarly via call_mutations_pileup.py
#to use from command line: 
#python pileup_parser.py mpileup_file > output_file

inFile = open(sys.argv[1],'r') #reads from stdin - fed by the call_mutations_pileup.py script

print 'CHROM\tpos\tref\tdepth\tA\tG\tC\tT\tN\tdel_count\tins_count\tdeleted\tinserted\tambiguous\tCall\tMutant' #header
for line in inFile: #parse lines in the input pileup file - each line indicates position in the reference sequence
	data = line.strip().split('\t') #split the input file by tabs, into a list
	while True: #loop through all characters in the line; breaks either after its completion or at the beginning if it fails certain checks pertaining to low coverage regions (no coverage or only 1 *)
	    if data[3]=='0': 
		break # break if the reported depth is 0 at this position
	    data_set=set(data[4].upper()) #a unique list of the bases
	    if len(data_set)<2 and '*' in data_set:
		break #break if there's only 1 unique base, and it's a * (placeholder for deletions)
	    pos = data[1] #assigning some variables based on the entries of the string. bases has what needs to be parsed
            bases = data[4].upper()
            ref = data[2].upper()
            chrom=data[0]
            types = {'A':0,'G':0,'C':0,'T':0,'N':0, '-':0,'+':0, '*':0, 'X':[], 'depth':0} # dictionary to store values for each all pileup variables (SNVs, occurrence of indels, ambiguous=X, depth) with specifics of indels stored in separate dictionaries delList and insList
	    delList={}
	    insList={}

            i = 0
	    bases = bases.replace('$', '') #$ signifies end of a read, but not necessary for parsing
            while i < len(bases): #length of bases is n, i ranges from 0:n-1
                    base = bases[i]
		    if base == '^': #new read character
                            i += 2 #^ is always followed another character. ascii of this character -33 = quality 
		    elif base == '-': # marks start of a deletion in pileup eg -2AG; start logging deletion specifics 
                            i += 1 #would move iterator to the '2' in the above example
			    types['-']+=1 #deletion count
			    indelLength = list() # number that indicates length of the indel (indels start with length specification in the pileup file)
                            while bases[i] not in ['A', 'C', 'G', 'T', ',', '.', '^', '-', '+', 'N']: #breaks when bases[i] equals various possible characters, most importantly the bases, and other characters which would signify the ent of the deltion. added N in case reference has ambiguous bases (e.g. CHO)
                                    temp = bases[i]
                                    indelLength.append(temp) #if the indel was something like -10AAAAAAAAAA, this appends the digits to a list i.e. ['1', '0']
				    i += 1
                            indelLength = ''.join(indelLength)
			    indelLength = int(indelLength) #in above example, concatenates digits and converts to int
			    delSeq = ''  # deleted sequence
			    for j in range(indelLength): #gathers the deleted sequence bases for the length of the indel
				    delSeq += bases[i]
				    i+=1
			    if delSeq in delList:   # delList is the dictionary of all unique deletions sequences and their counts (=occurrence at same position)  
				    delList[delSeq]+=1  # deletion is in the dictionary then add 1 to counts
			    else:  # if deletion not in dictionary, add it to the dictionary and set counts to 1
				    delList[delSeq]=1
                    elif base == '*': #the placeholder for deletions. e.g. if a line had a deletion of -1A, then the next line would have a *
                            types['*'] += 1 #currently counted, but not output
			    i+=1
	            elif base == '+': # marks start of a insertion in pileup file, then start logging insertion specifics 
                            i += 1
			    types['+'] +=1
			    indelLength = list()
			    while bases[i] not in ['A', 'C', 'G', 'T', ',', '.', '^', '-', '+', 'N']: #same as above for deletions. N accounts for insertions of ambiguous bases
			    	    temp = bases[i]
				    indelLength.append(temp)
				    i += 1
                            indelLength = ''.join(indelLength)
			    indelLength = int(indelLength)
			    insSeq = ''  # inserted sequence
                            for j in range(indelLength): #gathers the inserted sequence bases for the length of the indel
                                    insSeq += bases[i]
				    i+=1
			    if insSeq in insList: # insList is the dictionary of all unique insertion sequences and their counts (=occurrence at same position)
                                    insList[insSeq]+=1  # insertion is in the dictionary then add 1 to counts
                            else:  # if insertion is not in dictionary, add it to the dictionary and set counts to 1
                                    insList[insSeq]=1
                    elif base == '.' or base == ',': #ref base on forward or referse strand respectively
                            types[ref] += 1
			    types['depth'] +=1
			    i+=1
           	    else:   # if base is not reference or an indel 
                            if base in types: # check if base is an SNV or N
                                    types[base] += 1
				    types['depth'] +=1
				    i+=1
                            else: # else if base is no indel, SNV, N - it's specified as unidentified character which will be ignored and counted under ambiguous
                                    types['X'].append(base) # unidentified characters
				    i+=1
			    

            insOut = '' #"empty" character for fields. could also be ''
            if len(insList) > 0:
		    insOut = "; ".join(":".join((str(k),str(v))) for k,v in insList.items())
	    delOut = ''
	    if len(delList) > 0:
	            delOut = "; ".join(":".join((str(k),str(v))) for k,v in delList.items())
            amb = ''
            if len(types['X']) > 0:
                    amb = sorted(set(types['X']))  # make unique list of all unidentified characters encountered 
		    amb = ','.join(amb)

	    call=100*float(types[ref])/float(types['depth']) #calculating the rate of reference bases in read depth (counting only reference,SNV,N bases towards depth)
	    call = str(round(call, 4)) #rounding  to 4 decimal places
	    mutant = 100* (1-(float(types[ref]))/float(types['depth'])) #calculate mutation rate
	    mutant = str(round(mutant, 4)) # rounding
            out = [chrom, pos, ref, types['depth'],types['A'],types['G'],types['C'],types['T'],types['N'],types['-'],types['+'],delOut, insOut, amb, call, mutant] #big list
            print '\t'.join([str(x) for x in out]) #prints to stdout
	    # all pileup characters on line processed - breaking out of the line while loop
	    break
