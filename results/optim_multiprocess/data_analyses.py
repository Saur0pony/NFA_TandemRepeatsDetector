import matplotlib.pyplot as plt

chrom = []
length = []

times_mult1 = []
times_mult2 = []
times_mult3 = []

times_simp = []
times_simp1 = []
times_simp2 = []
times_simp3 = []

times_cpu2_1 = []
times_cpu2_2 = []
times_cpu2_3 = []

times_cpu6_1 = []
times_cpu6_2 = []
times_cpu6_3 = []


with open("comparaisonspacer.txt","r") as data:
    for line in data:
        values = line.split()
        chrom.append(values[0])
        length.append(int(values[1]))
        times_mult1.append(int(values[2]))
        times_mult2.append(int(values[3]))
        '''
        times_mult3.append(int(values[4]))
        times_simp1.append(int(values[5]))
        
        times_simp2.append(int(values[6]))
        times_simp3.append(int(values[7]))
        times_cpu2_1.append(int(values[8]))
        times_cpu2_2.append(int(values[9]))
        times_cpu2_3.append(int(values[10]))
        times_cpu6_1.append(int(values[11]))
        times_cpu6_2.append(int(values[12]))
        times_cpu6_3.append(int(values[13]))


data = [times_mult1, times_mult2, times_mult3]
times_mult = [round(sum(e)/len(e)) for e in zip(*data)]

data = [times_simp1, times_simp2, times_simp3]
times_simp = [round(sum(e)/len(e)) for e in zip(*data)]

data = [times_cpu2_1, times_cpu2_2, times_cpu2_3]
times_cpu2 = [round(sum(e)/len(e)) for e in zip(*data)]

data = [times_cpu6_1, times_cpu6_2, times_cpu6_3]
times_cpu6 = [round(sum(e)/len(e)) for e in zip(*data)]
'''
dicolen_mult = {}
dicolen_simp = {}
dicolen_cpu2 = {}
dicolen_cpu6 = {}


for i in range(len(length)):
    dicolen_mult[length[i]] = times_mult1[i]
    #dicolen_simp[length[i]] = times_simp1[i]
    dicolen_cpu2[length[i]] = times_mult2[i]
    #dicolen_cpu6[length[i]] = times_mult3[i]


true_length = []
true_times_mult = []
true_times_simp = []
true_times_cpu2 = []
true_times_cpu6 = []
for key in sorted(dicolen_mult.keys()):
    true_length.append(key)
    true_times_mult.append(dicolen_mult[key])
    #true_times_simp.append(dicolen_simp[key])
    true_times_cpu2.append(dicolen_cpu2[key])
    #true_times_cpu6.append(dicolen_cpu6[key])
    

#plt.plot(true_length, true_times_simp, label = '6CPUs - cut 5 Mb')
plt.plot(true_length, true_times_mult, label = '6CPUs - avec spacer')
plt.plot(true_length, true_times_cpu2, label = '6CPUs - sans spacer')
#plt.plot(true_length, true_times_cpu6, label = '6CPUs - cut 100 Mb')


plt.grid()
plt.xlabel("length (x10^8 b)")
plt.ylabel("time (s)")
plt.title("Execution time of pattern matching NFA in function of sequence length")
plt.legend()
plt.savefig("comp_sans_spacers.png")
