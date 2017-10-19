import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind_from_stats as ttest

def get_stats(data):
    std = np.std(data)
    mean = np.mean(data)
    return mean, std


# VUL IN DEZE TWEE EXCELBESTANDEN DE WAARDEN IN UIT DE TIEN RUNS PER ALGORITME #
# metrics: fitness increase per iteration and final fitness value
data = pd.read_excel('results_bench.xlsx')
data2 = pd.read_excel('results.xlsx')

# algo 1 (benchmark)
increase = np.asarray(data['fit_inc'])
fitness = np.asarray(data['fitness'])

# algo 2
increase2 = np.asarray(data2['fit_inc'])
fitness2 = np.asarray(data2['fitness'])

# algo 1 (benchmark)
m1, s1 = get_stats(increase[0:9])
m2, s2 = get_stats(increase[10:19])
m3, s3 = get_stats(increase[20:29]) 

M1, S1 = get_stats(fitness[0:9])
M2, S2 = get_stats(fitness[10:19])
M3, S3 = get_stats(fitness[20:29]) 

# algo 2 (benchmark)
m12, s12 = get_stats(increase2[0:9])
m22, s22 = get_stats(increase2[10:19])
m32, s32 = get_stats(increase2[20:29]) 

M12, S12 = get_stats(fitness2[0:9])
M22, S22 = get_stats(fitness2[10:19])
M32, S32 = get_stats(fitness2[20:29]) 

x = [.8,1.8,2.8] 
x2 = [1.2,2.2,3.2]  
xlabels = [1,2,3]  
labels = ['Bent cigar', 'Schaffers', 'Katsuura']    

# fitness increase
plt.subplot(121)
plt.errorbar(x, [m1,m2,m3], yerr=[s1,s2,s3], fmt='sr', ecolor='r', capsize=3)
plt.errorbar(x2, [m12,m22,m32], yerr=[s12,s22,s32], fmt='sb', ecolor='b', capsize=3)
ax = plt.gca()
ax.set_xticks(xlabels)
ax.set_xticklabels(labels)
# ax.annotate('local max', xy=(1, m1+s1), xytext=(1, m1+s1+1))
plt.ylim([0,3])
plt.xlim([0,4])
plt.grid()
plt.ylabel('fitness increase per iteration')

# total fitness
plt.subplot(122)
plt.errorbar(x, [M1,M2,M3], yerr=[S1,S2,S3], fmt='sr', ecolor='r', capsize=3)
plt.errorbar(x2, [M12,M22,M32], yerr=[S12,S22,S32], fmt='sb', ecolor='b', capsize=3)
ax = plt.gca()
ax.set_xticks(xlabels)
ax.set_xticklabels(labels)
plt.ylim([0,10])
plt.xlim([0,4])
plt.legend(['Benchmark (hillclimber)', 'Plant propagation algorithm'])
plt.grid()
plt.ylabel('Attained fitness')

plt.suptitle('Performance metrics for three evaluation functions')
plt.show()


## statistics: check significance of each test against the benchmark
_,p1 = ttest(m1, s1, 10, m12, s12, 10)
_,p2 = ttest(m2, s2, 10, m22, s22, 10)
_,p3 = ttest(m3, s3, 10, m32, s32, 10)

_,p4 = ttest(M1, S1, 10, M12, S12, 10)
_,p5 = ttest(M2, S2, 10, M22, S22, 10)
_,p6 = ttest(M3, S3, 10, M32, S32, 10)

print p1, p2, p3, p4, p5, p6