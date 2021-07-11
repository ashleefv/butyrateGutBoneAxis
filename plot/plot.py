import numpy as np
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import random,pickle

result = np.array(pickle.load(open('result.p', 'rb')))
result1 = np.array(pickle.load(open('result1.p', 'rb')))
result2 = np.array(pickle.load(open('result2.p', 'rb')))
result3 = np.array(pickle.load(open('resulth.p', 'rb')))
T = np.array(pickle.load(open('time.p', 'rb')))


plt.rcParams.update({'font.size': 25})
plt.xticks(np.arange(0, 30, step=7))


#plt.plot(T, result3[:,5]*100/(result3[:,2] + result3[:,5]) ,T, result[:,5]*100/(result[:,2] + result[:,5]))
#plt.plot(T, result3[:,6]*100/(result3[:,3] + result3[:,6]) ,T, result[:,6]*100/(result[:,3] + result[:,6]))
#plt.plot(T, (result3[:,7]/(result3[:,4] + result3[:,7]))*100 ,T, (result[:,7]/(result[:,4] + result[:,7]))*100)
plt.plot(T, result[:,0], T, result[:,1], linewidth=3)


#plt.fill_between(T, result1[:,5]*100/(result1[:,2] + result1[:,5]), result2[:,5]*100/(result2[:,2] + result2[:,5]), color = 'k', alpha = 0.1)
#plt.fill_between(T, result1[:,6]*100/(result1[:,3] + result1[:,6]), result2[:,6]*100/(result2[:,3] + result2[:,6]), color = 'k', alpha = 0.1)
#plt.fill_between(T, (result1[:,7]/(result1[:,4] + result1[:,7]))*100, (result2[:,7]/(result2[:,4] + result2[:,7]))*100, color = 'k', alpha = 0.1)
plt.fill_between(T, result1[:,0], result2[:,0], color = 'k', alpha = 0.1)
plt.fill_between(T, result1[:,1], result2[:,1], color = 'k', alpha = 0.1)



#plt.plot(T, result[:,0]/volLymphT, T, result[:,1]/volLymphT,T, result[:,2]/volBlood, T, result[:,3]/volBone, T, result[:,4]/volBone)
#plt.plot(T, result[:,0], T, result[:,1],T, result[:,2], T, result[:,3], T, result[:,4])
#plt.yticks(np.arange(41, 48, step=1))
#plt.plot(T, (result[:,1]/(result[:,7] + result[:,1]))*100,T, (result1[:,1]/(result1[:,7] + result1[:,1]))*100, linewidth=3)

#plt.plot(T, (result[:,2]/(result[:,7] + result[:,2]))*100,T, (result1[:,2]/(result1[:,7] + result1[:,2]))*100, linewidth=3)
#plt.plot(T, (result[:,3]/(result[:,7] + result[:,3]))*100,T, (result1[:,3]/(result1[:,7] + result1[:,3]))*100, linewidth=3)
#plt.plot(T, result[:,3]/(result[:,8] + result[:,3]))
#plt.legend(['Control', 'Control + butyrate'], loc='upper left')
plt.legend(['butyrate in intestine', 'butyrate in blood'], loc='lower right')
plt.xlabel('TIME (days)')
plt.ylabel('$\Delta$ Butyrate \u03bcM')
#plt.ylabel('Bone Tregs (%)')
#plt.errorbar([28, 28],[42, 48], yerr=[1,2], fmt='o',color='r', elinewidth=3, markersize=10, capsize=6, capthick=3, barsabove= False)
#plt.errorbar([28, 28],[13, 16], yerr=[0.6,0.4], fmt='o',color='r', elinewidth=3, markersize=10, capsize=6, capthick=3, barsabove= False)
#plt.plot(T, result[:,3])
#plt.ylim([42,47.5])



plt.show()