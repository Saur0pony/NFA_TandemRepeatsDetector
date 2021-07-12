import getopt
import sys, os
import time
import multiprocessing as mp
from multiprocessing import Process, Pool, Queue, cpu_count
from threading import Thread


version = "NFA pattern matching v5 - no spacer vector - Multiprocessing"

"""Version note:
Modification of D update to take off "0" spacers
"""

def make_gapped_pattern(pattern, gap):
    """Function to make the gapped pattern for tandem repeat analyse

    Args:
        pattern (str): Pattern sequence
        gap (int): length of the gap between 2 pattern

    Returns:
        pattern_inti (str): Return the modified pattern
    """
    pattern_new = pattern
    for _ in range(gap):
        pattern_new += "N"
    pattern_new += pattern
    return pattern_new


def preprocessing(pattern, m, k):
    """Function use to create initial state of variables which are bit arrays represented by correspondant integer

    Args:
        pattern (str): pattern for the pattern matching
        m (int): length of the pattern
        k (int): number of maximum errors allowed

    Returns:
        T, S, M1, M2, Din, G (int): return all the constant we need + the initial state of diagonal binary vector D in int
    """
    t = {}
    T = {}
    S = {}
    alphabet = ["A", "C", "G", "T", "N", "$", "R", "Y", "M", "W", "S", "K", "B", "D", "H", "V", "O", "I", "F"]

    for char in alphabet:
        t[char] = 0
        for i in range(m):
            if char != pattern[i] and pattern[i] != "N":
                t[char] += 2**i
        #print(f't({char}) = {t[char]}')

        if char == "N":
            t[char] = 0

        T[char] = 0
        for i in range(m-k):
            # get k+1 last bits of t(char)
            T_inter = ((t[char] >> i) & (2**(k+1)-1))
            T[char] += T_inter << (m-k-1-i)*(k+1)

        #print(f'T({char}) = {T[char]} ')
        S[char] = char in pattern[:k + 1]

    # constant mask 1 : sub / del (for diagonal transition)
    M1 = 0
    for i in range(0, (m-k)*(k+1), k+1): 
        M1 += int(2**(i))
    print(f'M1 = {M1}')

    # constant mask 2 : insert (for vertical transition)
    M2 = 2**(k+1) - 2
    for i in range(0, (m-k)*(k+1), k+1):
        M2 += int(2**(i))
    print(f'M2 = {M2}')

    # diagonals bit array set to 1 for all states of all diags = initial state
    Din = (2**((k+1)*(m-k)))-1
    #print(f'Din = {Din}')

    # constant to calculate z:
    Mz = 0
    One = ( 1<<k ) - 1
    for i in range(m-k):
        Mz = (Mz << (k+1)) + One

    G = 1 << k

    Gbis = (2**(k+1))-1

    return T, S, M1, M2, Mz, Din, G, Gbis


def update_d(D, Tlettre, k, M1, M2, Mz):
    """Base function of NFA algo to update the D vector

    Args:
        D (int): The D vector that contain all state by diagonals
        Tlettre (int): Match mask for the actual letter at the actual position
        k (int): max errors allowed number
        M1 (int): Mask 1 (0**k 1)**m-k 
        M2 (int): Mask 2 (0**k 1)**m-k-1 1**k+1
    """
    x = (D >> (k + 1)) | Tlettre
    z = x & Mz
    
    """
    return(
        (((D << 1) | M1) & ((D << (k+2)) | M2)) 
        & (((z + M1) ^ z) & x))
    """

    '''
    return(
        ((D << 1) | M1)  # substitution
        & ((D << (k+3)) | M2)  # insertion
        & ((x + M1) ^ x) >> 1  # match and deletion
        & Din)
    '''
    #print(f'x = {x}')

    return(
        (((D << 1) & ((z + M1) ^ z)) | M1)
        & ((D << (k+2)) | M2) & x
    )


def get_nb_error(D, Gbis):
    """Function to get the error number of the last diagonal match

    Args:
        D (int): The diagonals vector
        Gbis ([type]): mask to get last diagonal in the vector

    Returns:
        nb_error (int): Calculated number of errors at the match position
    """
    nb_error = 0
    Dk = D & Gbis
    while Dk:
        Dk = Dk & (Dk-1)
        nb_error += 1
    return nb_error


def pattern_search(pattern, sequence, k, gap, cutpat, ID):
    """Function use to get list of occurence of pattern in a sequence. It use binary optimisation with bitarray that correspond
    to diagonale states. Those are called vector and are represented with an integer.

    Args:
        pattern (str): The wanted pattern sequence 
        sequence (str): The sequence where we want to search the pattern
        k (int): max allowed error in the pattern / word distance

    Returns:
        list_occur(list): list with all positions where the pattern was found
    """
    starting_time = time.time()
    n = len(sequence)
    m = len(pattern)

    '''Preprocessing'''
    #for the initial pattern
    T_init, S_init, M1_init, M2_init, Mz_init, Din_init, G, Gbis = preprocessing(pattern, m, k)
    initial_pattern = pattern
    m_init = m

    #for the (N)n-pattern
    cut_pattern = pattern[-int(m-((m-gap-k)/2)):]
    m_cut = len(cut_pattern)
    T_cut, S_cut, M1_cut, M2_cut, Mz_cut, Din_cut, G, Gbis = preprocessing(cut_pattern, m_cut, k)

    T, S, M1, M2, Mz, Din = T_init, S_init, M1_init, M2_init, Mz_init, Din_init
    tandem_repeat = False

    '''searching for pattern occurences'''
    list_occur = []
    sub_list_occur = []
    D = Din
    i = 0
    while i < n:
        if (S[sequence[i]] or tandem_repeat) and sequence[i] != "N":
            while i < n:
                #print(f'==== Pos {i} : {sequence[i]} ====')

                D = update_d(D, T[sequence[i]], k, M1, M2, Mz)
                #print(D)

                if D & G == 0:
                    #error count
                    nb_error = get_nb_error(D, Gbis)
                    #print(nb_error)
                    #true i identification and save position
                    i -= nb_error
                    sub_list_occur.append(i+ID)
                    #print(i+ID)

                    #keep  only gap_pattern if it was the first occurence to elong the TR
                    if cutpat and pattern == initial_pattern:
                        pattern = cut_pattern
                        m = m_cut
                        T, S, M1, M2, Mz, Din = T_cut, S_cut, M1_cut, M2_cut, Mz_cut, Din_cut
                        tandem_repeat = True
                    
                    D = Din  # reset to search a new pattern
                      
                #check if we are in TR and if the actual position is far from last occurence saved position
                elif tandem_repeat and i+ID - sub_list_occur[-1] > m+k :
                    pattern = initial_pattern
                    m = m_init
                    T, S, M1, M2, Mz, Din = T_init, S_init, M1_init, M2_init, Mz_init, Din_init
                    tandem_repeat = False
                    D = Din
                    list_occur.append(sub_list_occur)
                    sub_list_occur = []

                if D == Din:
                    # if we come back to initial vector because of no more matching, break
                    break
                i += 1
        i += 1
    if sub_list_occur != []:
        list_occur.append(sub_list_occur)

    print(starting_time - time.time())

    print(list_occur)

    return list_occur


def worker(in_q, out_q, pattern, k ,gap, cutpat):
    while True:
        try:
            sequence, pos_init = in_q.get(True)
            out_q.put(pattern_search(pattern, sequence, k, gap, cutpat, pos_init))
        except TypeError:
            in_q.put(None)
            return

def outworker(out_queue, out_list):
    while True:
        line = out_queue.get(True)
        if line is None:
            break
        for e in line:
            print(e)
        out_list += line

def multi_process_cut(pattern, sequence, k, gap, cutpat, step):
    
    nb_core = os.cpu_count()
    print(f'nb of CPU : {nb_core}')

    out_list = []
    #($)k to add at the end of each sequence for PM
    addit = ""
    for _ in range(k):
        addit += "$"

    #cut sequence into pieces
    sequences = []
    start = 0
    end = int(step/500) + step
    while end < len(sequence)+k:
        sequences.append(sequence[start:end] + addit)
        start += step
        end += step
    sequences.append(sequence[start:])
    

    in_q = Queue(maxsize=nb_core)
    out_q = Queue(maxsize=200)

    writer = Thread(target=outworker, args=(out_q, out_list))
    writer.start()

    processes = Pool(processes=nb_core, initializer=worker, initargs=(in_q, out_q, pattern, k, gap, cutpat))

    pos_init = 0
    for seq in sequences:
        in_q.put((seq, pos_init))
        pos_init += step

    
    in_q.put(None)
    processes.close()
    processes.join()
    out_q.put(None)
    writer.join()

    #print(out_list)
    return(out_list)


def get_arguments():
    """function to get arguments from execute command 
    --file : fasta file that contain sequences to analyse
    --out : name of output file with results
    --pattern : the pattern to search in the sequence
    --k : value of max error allowed

    Returns:
        file_name (str) : name of the fasta file that contain sequences
        output_name (str) : name of the otput file
        pattern (str) : the pattern to search in sequences
        k (int) : maximum error allowed
        gap (int) : gap length between 2 motifs
    """
    file_name = None
    output_name = None
    pattern = None
    k = 2  # define k standard value to 2
    gap = 0
    cutpat = True

    try:
        options, _ = getopt.getopt(
            sys.argv[1:], "", ["file=", "out=", "pattern=", "k=", "gap=", "cutpat", "help"])
    except getopt.GetoptError as err:
        print(f'{err}. Please use agument --help to see how to execute the pattern_matching_NFA script')
        sys.exit(2)

    '''getting values from arguments'''
    for opt, arg in options:
        if opt in "--file":
            file_name = arg
        elif opt in "--out":
            output_name = arg
        elif opt in "--pattern":
            pattern = arg
        elif opt in "--k":
            k = int(arg)
        elif opt in "--gap":
            gap = int(arg)
        elif opt in "--cutpat":
            cutpat = False

    return file_name, output_name, pattern, k, gap, cutpat


def main():
    """main function"""
    start_time = time.time()

    #get arguments
    file_name, output_name, pattern, k, gap, cutpat = get_arguments()


    if gap != None:
        pattern = make_gapped_pattern(pattern, gap)

    extend = ""
    for _ in range(k):
        pattern += "N"
        extend += "$"

    read_seq = False
    id_seq = ""
    list_pos = []
    with open(output_name, "w") as output:
        output.write(
            f'>Pattern matching by using NFA algorithme : {version}\n>Sequence file name : {file_name}\n>Pattern : {pattern}\n>k = {k}\n')

        with open(file_name, "r") as sequence:
            for line in sequence:
                if not read_seq and line.startswith(">"):
                    id_seq = line.rstrip()
                    read_seq = True
                elif read_seq:
                    #list_pos = pattern_search(pattern, line.rstrip() + extend, k, gap, cutpat)
                    list_pos = multi_process_cut(pattern, line.rstrip() + extend, k, gap, cutpat, 10000000)
                    output.write(f'>{id_seq}\n>Length:\t{len(line.rstrip())}\n>Bloc_found:\t{len(list_pos)}\n>Execution_time:\t{round(time.time() - start_time,2)}\n')
                    output.write(f'>start\tend\tcopies\tlength\tpositions\n')
                    nb_copy = 0
                    for tandem in list_pos:
                        output.write(f'{tandem[0] - (len(pattern)-k)}\t{tandem[-1]}\t{len(tandem)+1}\t{tandem[-1] - tandem[0] + len(pattern)-k}\t{tandem}\n')
                        nb_copy += len(tandem)+1
                    if list_pos != []:
                        output.write(f'Global :\n{list_pos[0][0] - (len(pattern)-k)}\t{list_pos[-1][-1]}\t{nb_copy}\t{list_pos[-1][-1] - (list_pos[0][0] - (len(pattern)-k))}\n')

                    read_seq = False


if __name__ == "__main__":

    
    print("Starting !!!")
    start_time = time.time()

    main()

    total_time = round(time.time() - start_time,2)
    m_total_time = int(total_time // 60)
    s_total_time = total_time % 60

    print(f'>Execution Time : 00:{m_total_time:02d}:{int(s_total_time//1):02d}.{int(round(s_total_time%1,2)*100):02d}')
    print("End of process !!!") 

    

    """ DEBUG SECTION : Some Test for the program """

    # execution test

    '''pattern = "ATTAAGGG"
    kk = 2
    gagap = 2
    cucutpat = True
    sequence = "CCCCATTAGGGCCATTAAGGGCCATTAAGGGCCATTAAGGGCCATTAAGGG"

    strgap = ""
    for _ in range(gagap):
        strgap += 'N'

    pattern += strgap + pattern
    extand =""
    for _ in range(kk):
        pattern += "N"
        extand += "$"

    positions = pattern_search(pattern, sequence + extand, kk, gagap, cucutpat)
    print(positions) 
    '''
    '''
    multi_process_cut("ATTAAGGGNN",   "ACTGTTGTACACATGTTGTCCCCACACACTGTTGTACACACAATTAAGGGATTAAGGGATTAAGGGCACACTGTATGTCACACCATTAAGGGATTAAGGGCCCTGTTGTTGTCACACACACACTGTTGTGCACACACACACTTTGTGTGGGTGTGTTGTTGTTGTACATGTGTGTGTGTTGTGTTGTTGTTACACTGATTAAGGGATTAAGGGATTAAGGGATTAAGGGACTACTGTGTGAC$$", 2, 0, True, 200)
    
    #print(the_liste)
    queue = Queue()
    #the_other_list = pattern_search("TGTTGTN", "ACTGTTGTACACATGTTGTCCCCACACACTGTTGTACACACACACACTGTATGTCACACCCCCTGTTGTTGTCACACACACACTGTTGTGCACACACACACTTTGTGTGGGTGTGTTGTTGTTGTACATGTGTGTGTGTTGTGTTGTTGTTACACTGACTACTGTGTGAC$", 1, 0, True, queue, 0)
    #print(the_other_list)

    pattern_search("ATTAAGGGNN",      "ACTGTTGTACACATGTTGTCCCCACACACTGTTGTACACACAATTAAGGGATTAAGGGATTAAGGGCACACTGTATGTCACACCATTAAGGGATTAAGGGCCCTGTTGTTGTCACACACACACTGTTGTGCACACACACACTTTGTGTGGGTGTGTTGTTGTTGTACATGTGTGTGTGTTGTGTTGTTGTTACACTGATTAAGGGATTAAGGGATTAAGGGATTAAGGGACTACTGTGTGAC$$", 2, 0, True, 0)
    '''