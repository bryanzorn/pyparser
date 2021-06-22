import re

SampleDict = {}
with open("samplelist.txt") as samples:

    for line in samples:
        (key, val) = line.split()
        SampleDict[key] = val

bcodes = [line.rstrip('\n') for line in open('barcodes.txt')]

for i in range(len(bcodes)):

	curcode = str(bcodes[i].strip('\r'))

	with open('fastqjoin.join') as data:

		outfile = open("demultiplexed_fastqjoin.join_seqs.fna", 'a+')
		seqNumb = 0

		for line in data:
			seqName = line.strip('\n')
			seqInfo = seqName[0:-8]

			illumbc = seqName[-8:]
			demuxbc = str(curcode+illumbc)

			seqData = data.next().strip('\n')
                        inlineF = 4+int(len(curcode))+4
                        inlineR = 2+(8-int(len(curcode)))+3

			plusSep = data.next().strip('\n')
			qScores = data.next().strip('\n')

			if re.match("...." + curcode + "...." + "TCACTCCTACGGGAGG", seqData):
                                for k,v in SampleDict.items():
                                    if demuxbc == v:
                                        SampleID = k
				seqNumb += 1
				outfile.writelines([">", str(SampleID), "_", str(seqNumb), \
                                " ", seqInfo, " new_bc=", demuxbc, " bc_diffs=0", "\n", \
                                seqData[inlineF:-inlineR], "\n"])
		outfile.close()
