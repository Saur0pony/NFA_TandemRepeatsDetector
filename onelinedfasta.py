import sys


with open("GCF_000001405.39_GRCh38.p13_genomicmodifiiii.fna","w") as output:
    with open("GCF_000001405.39_GRCh38.p13_genomic.fna","r") as file:
        sequence = ""
        start = False
        needseq = False
        for line in file:
            if ">" in line:
                if ">NC_" in line:
                    
                    if start and needseq:
                        sequence += "\n"
                        output.write(sequence)
                        #print(f'{sequence[0:9]} [...] {sequence[-10:]}')
                        sequence = ""
                    
                    needseq = True

                    output.write(line)
                    print(line.rstrip())
                    
                else:
                    if start and needseq:
                        sequence += "\n"
                        output.write(sequence)
                        #print(f'{sequence[0:9]} [...] {sequence[-10:]}')
                        sequence = ""
                    needseq = False
                       
                start = True

            elif needseq:
                sequence += line.rstrip().upper()
        if needseq:     
            output.write(sequence)
