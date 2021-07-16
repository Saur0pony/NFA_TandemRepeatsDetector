import sys, os, getopt


def get_arguments():
    """function to get arguments from execute command 
    --file : fasta file that contain sequences to transform
    --out : name of output file with results

    Returns:
        file_name (str) : name of the fasta file that contain sequences
        output_name (str) : name of the otput file
    """

    file_name = None
    output_name = None

    try:
        options, _ = getopt.getopt(
            sys.argv[1:], "", ["file=", "out=", "help"])
    except getopt.GetoptError as err:
        print(f'{err}. Please use agument --help to see how to execute the Maxrep script')
        sys.exit(2)

    '''getting values from arguments'''
    for opt, arg in options:
        if opt in "--file":
            file_name = arg
        elif opt in "--out":
            output_name = arg

    return file_name, output_name


if __name__ == "__main__":

    file_name, output_name = get_arguments()

    with open(output_name,"w") as output:
        with open(file_name,"r") as file:
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
