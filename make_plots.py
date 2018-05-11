import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


names = np.chararray( 24,unicode=True, itemsize=15)
percent = np.zeros(24)

f = open('output_c3.txt')
i = 0
for line in f.readlines():
    # print(line)
    templine = line.split('\t')
    # print(templine)
    temp = templine[0]
    temp = temp.strip('.TIF').split('_')
    names[i] = temp[1]
    percent[i] = float(templine[3])
    print(names[i],percent[i])
    i+=1

PTX20 = percent[[1,13]]
PTX10 = percent[[3,15]]
TGFB2 = percent[[5,17]]
TGFB1 = percent[[7,19]]
untreated = percent[[9,21]]
stopper = percent[[11,23]]

PTX20 = percent[[0,12]]
PTX10 = percent[[2,14]]
TGFB2 = percent[[4,16]]
TGFB1 = percent[[6,18]]
untreated = percent[[8,20]]
stopper = percent[[10,22]]

x = np.array([1,2,3,4,5,6])
width = 0.8
means = [np.average(untreated),np.average(stopper),np.average(TGFB1),np.average(TGFB2),np.average(PTX10),np.average(PTX20)]
var =  [np.var(untreated),np.var(stopper),np.var(TGFB1),np.var(TGFB2),np.var(PTX10),np.var(PTX20)]
sem = [stats.sem(untreated),stats.sem(stopper),stats.sem(TGFB1),stats.sem(TGFB2),stats.sem(PTX10),stats.sem(PTX20)]

for i in [0,1]:
    plt.errorbar(x=1+width*0.125+width*(i*0.25),y=untreated[i],fmt='ko')
    plt.errorbar(x=2+width*0.125+width*(i*0.25),y=stopper[i],fmt='ko')
    plt.errorbar(x=3+width*0.125+width*(i*0.25),y=TGFB1[i],fmt='ko')
    plt.errorbar(x=4+width*0.125+width*(i*0.25),y=TGFB2[i],fmt='ko')
    plt.errorbar(x=5+width*0.125+width*(i*0.25),y=PTX10[i],fmt='ko')
    plt.errorbar(x=6+width*0.125+width*(i*0.25),y=PTX20[i],fmt='ko')

# plt.bar(left=x,height=means,yerr=sem,width=width)
plt.bar(left=x,height=means,width=width)
ax = plt.gca()
ax.set_xticks(x)
ax.set_xticklabels(('Untreated', 'Stopper', 'TGFB1', 'TGFB2','PTX10','PTX20'))
plt.xlim(xmin=0.2,xmax=6+width)

vals = ax.get_yticks()
ax.set_yticklabels(['{:4.0f}%'.format(x) for x in vals])


plt.ylabel('Migration region covered by cells')
plt.xlabel('Experimental condition')
plt.savefig('comparison.png')

#0=untreated
#1=stopper
#2=TGFB1
#3=TGFB2
#4=PTX10
#5=PTX20
n=2 #two replicates
SE_TGFB1 = np.sqrt(var[0]/n+var[2]/n)
DF_TGFB1 = np.round((var[0]/n+var[2]/n)**2 / ( (var[0]/n)**2/3 + (var[2]/n)**2/3 ))
t_TGFB1 = (means[0]-means[2])/SE_TGFB1
p_TGFB1 = stats.t.cdf(x=t_TGFB1,df=DF_TGFB1)
print('TGFB1')
print(SE_TGFB1,DF_TGFB1,t_TGFB1,p_TGFB1)

SE_TGFB2 = np.sqrt(var[0]/n+var[3]/n)
DF_TGFB2 = np.round((var[0]/n+var[3]/n)**2 / ( (var[0]/n)**2/3 + (var[3]/n)**2/3 ))
t_TGFB2 = (means[0]-means[3])/SE_TGFB2
p_TGFB2 = stats.t.cdf(x=t_TGFB2,df=DF_TGFB2)
print('TGFB2')
print(SE_TGFB2,DF_TGFB2,t_TGFB2,p_TGFB2)

SE_PTX10 = np.sqrt(var[0]/n+var[4]/n)
DF_PTX10 = np.round((var[0]/n+var[4]/n)**2 / ( (var[0]/n)**2/3 + (var[4]/n)**2/3 ))
t_PTX10 = (means[4]-means[0])/SE_PTX10
p_PTX10 = stats.t.cdf(x=t_PTX10,df=DF_PTX10)
print('PTX10')
print(SE_PTX10,DF_PTX10,t_PTX10,p_PTX10)

SE_PTX20 = np.sqrt(var[0]/n+var[5]/n)
DF_PTX20 = np.round((var[0]/n+var[5]/n)**2 / ( (var[0]/n)**2/3 + (var[5]/n)**2/3 ))
t_PTX20 = (means[5]-means[0])/SE_PTX20
p_PTX20 = stats.t.cdf(x=t_PTX20,df=DF_PTX20)
print('PTX20')
print(SE_PTX20,DF_PTX20,t_PTX20,p_PTX20)

print(means)
print(sem)

