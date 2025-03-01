#!/usr/bin/python3 -u
import math
import sys,gzip,os
from os.path import basename
from statistics import variance

from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.SeqUtils import gc_fraction #GC #gc_fraction
import time

#start
start_time = time.time()

arguments = sys.argv[1:]

print("File\tA_percent\tT_percent\tC_percent\tG_percent\tGC_percent\tAT/GC_ratio\tNucleScore\tATG\tTGA\tTAG\tTAA\tGenome_size\tSpecies")

for argv in arguments:
    #print("File: "+argv+"\n\n")
    file = argv
    globalSeq = ''
    fasta_sequences = ''

    if file.endswith(".gz"):
        with gzip.open(file, 'rt') as fastin:
            with open('tmp.fna', 'wt') as fastout:
                data = fastin.read()
                fastout.write(data)
                file = "tmp.fna"
                fastout.close()
        #fasta_sequences = list(SeqIO.parse(gzip.open(file), 'fasta'))
        #print("gz")
    
    fasta_sequences = list(SeqIO.parse(open(file), 'fasta'))

    description = fasta_sequences[0].description.split()
    species = "_".join(description[1:5])
    if "variant" in description:
        species = "_".join(description[1:5])
    else:
        species = "_".join(description[1:3])

    for seq_record in fasta_sequences:
        globalSeq += str(seq_record.seq)

    #gcpercent = gc_fraction(globalSeq) * 100
    
    my_dna = Seq(globalSeq)

    ade = my_dna.count("A")
    thy = my_dna.count("T")
    gua = my_dna.count("G")
    cyt = my_dna.count("C")
    n = my_dna.count("N")
    length = len(globalSeq)

    gcpercent = (gua+cyt)/(ade+thy+gua+cyt) * 100
 
    aPercent = (ade/length)*100
    tPercent = (thy/length)*100
    gPercent = (gua/length)*100
    cPercent = (cyt/length)*100
    nPercent = (n/length) * 100

    atgcRatio = (ade + thy) / (gua + cyt)
    percentList = (aPercent, tPercent, gPercent, cPercent, nPercent)
    variance_value = variance(percentList)
    nucleScore = math.log2((variance_value * gcpercent * atgcRatio ** 3) / math.sqrt(length))

    # act = my_dna.find("ACT")

    atg = my_dna.count('ATG')
    tga = my_dna.count('TGA')
    tag = my_dna.count('TAG')
    taa = my_dna.count('TAA')

    label = basename(argv)

    # Summary file
    print(label + "\t" + str(
            aPercent) + "\t" + str(tPercent) + "\t" + str(cPercent) + "\t" + str(gPercent) + "\t" + str(
            gcpercent) + "\t" + str(atgcRatio) + "\t" + str(nucleScore) + "\t" + str(atg) + "\t" + str(tga) + "\t" + str(tag) + "\t" + str(taa) + "\t" + str(length) + "\t" + species)

    if os.path.exists('./tmp.fna'):
        os.remove('./tmp.fna')
#end
end_time = time.time()
print(f"Total computation time: {end_time - start_time:.2f} seconds")
