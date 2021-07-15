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


with open("comparaison_spacer_k1-9.txt","r") as data:
    for line in data:
        values = line.split()
        chrom.append(int(values[0]))
        length.append(int(values[1]))
        times_mult1.append(90.34/float(values[2]))
        times_mult2.append(90.34/float(values[3]))
        times_mult3.append(248.96/float(values[4]))
        times_simp1.append(248.96/float(values[5]))
        '''
        
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
    

#plt.plot(true_length, true_times_simp, label = '6CPUs - cut 5 Mb')
plt.plot(chrom, times_mult1, label = 'Ch16 - v5 - sans spacer')
plt.plot(chrom, times_mult2, label = 'Ch16 - v4 - avec spacer')
plt.plot(chrom, times_mult3, label = 'Ch1 - v5 - sans spacer')
plt.plot(chrom, times_simp1, label = 'Ch1 - v4 - avec spacer')

#plt.plot(true_length, true_times_cpu6, label = '6CPUs - cut 100 Mb')


plt.grid()
plt.xlabel("k value")
plt.ylabel("Mb/sec")
plt.title("Execution time for Human Chromosome depending on k value")
plt.legend()
plt.savefig("comparaison_human_k1-10.png")
