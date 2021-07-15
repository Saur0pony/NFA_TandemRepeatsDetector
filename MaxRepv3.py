import sufarray
import time
import getopt, sys
#import bisect

from numpy import argsort, empty_like, arange


def bisect_left(a, x, lo=0, hi=None):
    """Return the index where to insert item x in list a, assuming a is sorted.
    The return value i is such that all e in a[:i] have e < x, and all e in
    a[i:] have e >= x.  So if x already appears in the list, a.insert(x) will
    insert just before the leftmost x already there.
    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.
    """

    if lo < 0:
        raise ValueError('lo must be non-negative')
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        # Use __lt__ to match the logic in list.sort() and in heapq
        if a[mid] < x: lo = mid+1
        else: hi = mid
    return lo

def find_lt(a, x):
    'Find rightmost value less than x'
    i = bisect_left(a, x)
    if i:
        return a[i-1] ,i-1
    raise ValueError

def get_sa(w):
    return sufarray.SufArray(w).get_array()

def get_inv(sa):
    inv = empty_like(sa)
    inv[sa] = arange(len(inv), dtype=inv.dtype)
    return inv

def get_lcp(sequence, sa, inv):
    n = len(sa)
    l = 0
    LCP = [0] * n
    for i in range(n-1):
        if inv[i] > 0:
            k = inv[i]
            j = sa[k-1]
            while sequence[i+l] == sequence[j+l]:
                l += 1
            LCP[inv[i]] = l
            if l > 0:
                l = l-1
    return LCP

def get_LCP(w, sa):
    n = len(w)
    lcp = []
    for i in range(len(sa) - 1):
        i0 = sa[i]
        i1 = sa[i+1]
        lcp.append(0)
        for j in range(n-max(i0,i1)):
            if w[i0 + j] == w[i1 + j]:
                lcp[-1] += 1
            else:
                break
    return lcp

def get_sa_table(sequence):
    print("Je calcul la SA")
    sa = get_sa(sequence)
    print("Je calcul la INVsa")
    inv = get_inv(sa)
    print("je calcul les LCP")
    lcp = get_LCP(sequence, sa)
    print("fini pour le tableau des Suffixes")
    return sa, inv, lcp

def find_maxrep(w, ml=5, output="output.out"):
    n = len(w)
    r, p, LCP = get_sa_table(w)
    S = [u for u in range(len(LCP)) if LCP[u] < ml]
    S.sort()
    S = [-1] + S + [n-1]
    I = argsort(LCP)

    count = 0

    with open(output,"w") as outfile:

        for t in range(min([t for t in range(len(I)) if LCP[I[t]] >= ml]), n-1):
            i   = I[t]
            value, value_index = find_lt(S,i)
            p_i = value + 1
            n_i = S[value_index + 1]
            S.insert(value_index + 1,i)
            if (p_i==0 or LCP[p_i-1]!=LCP[i]) and (n_i==n-1 or LCP[n_i]!=LCP[i]):                          # Maximal on right?
                if r[p_i]==0 or r[n_i]==0 or w[r[p_i]-1]!=w[r[n_i]-1] or p[r[n_i]-1]-p[r[p_i]-1]!=n_i-p_i: # Maximal on left?
                    count += 1
                    #print(f"J'en ai {count} !!!")
                    outfile.write(w[r[i]:r[i]+LCP[i]])
                    outfile.write("\n")

def cut_seq(seq, step):
    sequences = []
    start = 0
    end = int(step/500) + step
    while end < len(seq):
        sequences.append(seq[start:end] + "$")
        start += step
        end += step
    sequences.append(seq[start:] + "$")
    return sequences


def get_arguments():
    """function to get arguments from execute command 
    --file : fasta file that contain sequences to analyse
    --out : name of output file with results
    --ml : minimal length of patterns

    Returns:
        file_name (str) : name of the fasta file that contain sequences
        output_name (str) : name of the otput file
        ml (int) : value of minimal length allowed
    """

    file_name = None
    output_name = None
    ml = 20
    maxl = 100

    try:
        options, _ = getopt.getopt(
            sys.argv[1:], "", ["file=", "out=", "ml=", "maxl=", "help"])
    except getopt.GetoptError as err:
        print(f'{err}. Please use agument --help to see how to execute the Maxrep script')
        sys.exit(2)

    '''getting values from arguments'''
    for opt, arg in options:
        if opt in "--file":
            file_name = arg
        elif opt in "--out":
            output_name = arg
        elif opt in "--ml":
            ml = int(arg)
        elif opt in "--maxl":
            maxl = int(arg)

    return file_name, output_name, ml


if __name__ == "__main__":
    
    print("je commence")
    file_name, output_name, ml = get_arguments()

    print("j'ai bien pris les arguments")
    
    if file_name != None:
        print("il y a bien un fichier")
        with open(file_name, 'r') as file:
            for line in file:
                print("c'est une ligne")
                if ">" not in line:
                    seq = line
                    print("AAAAAaaaaaahhhh")
                    find_maxrep(seq + "$", ml, output_name)
                    
    else:
        find_maxrep("CCCATGAGTAGTGATGAGCGTAGTGAGTGACTGTCGCTGCTAGCTGGCATGCTAGTCAGGGGGTGTCGAGTCGATGCACTACAGAATATAGTGATGATGATGATGATTAAGGGATTAAGGGACGTGCAGCTGACTGACTGGATCGACTGATCGGATCGACTGAGTCAGTCGACTGCCCATTAAGGGATTAAGGGATTAAGGGCGTGCGTGTCGCCCCCGCGCGGGTTGTGTGTCGTCGTCGTCGGTCGTCGGTTGCTTGTTGCTTGTTGCTTGTTGCTT$", ml, output_name)

