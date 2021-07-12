import math
import numpy as np
import random,pickle
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

sc = 14  # scaling factor
cyclelength=sc
Nc=2

result_f1 = np.array(pickle.load(open('case1.p', 'rb')))     # no butyrate
result_f2 = np.array(pickle.load(open('case2.p', 'rb')))     # with butyrate + carley
result_f3 = np.array(pickle.load(open('case3.p', 'rb')))    # with butyrate + carley + osteoblast proliferation
result_f4 = np.array(pickle.load(open('case4.p', 'rb')))    # with butyrate + carley + butyrate increase pre-osteoblast to osteoblast differentiation
result_f5 = np.array(pickle.load(open('case5.p', 'rb')))    # with butyrate + carley + butyrate increase wnt10b dependent pre-osteoblast to osteoblast differentiation
result_f6 = np.array(pickle.load(open('case6.p', 'rb')))    # with butyrate + carley + butyrate inhibiting osteoclast by osteoclast death
result_f7 = np.array(pickle.load(open('case7.p', 'rb')))    # with butyrate + carley + butyrate inhibiting osteoclast differentiation
result_f8 = np.array(pickle.load(open('case8.p', 'rb')))    # with butyrate + carley + butyrate increase osteoblast mediated bone formation
result_f9 = np.array(pickle.load(open('case9.p', 'rb')))    # with butyrate + carley + butyrate reduce osteoclast mediated bone resorption
T1 = np.arange(0.0, Nc*cyclelength, 0.001)

#print(((np.trapz(y=result_f3[:,13], x=T1))-(np.trapz(y=result_f1[:,13], x=T1)))/(np.trapz(y=result_f1[:,13], x=T1)))
#print((np.trapz(y=result_f3[:,13], x=T1)), (np.trapz(y=result_f1[:,13], x=T1)))



# Osteoblast/Osteoclast area under curve
AUC1 = (np.trapz(y=result_f1[:,13], x=T1))/(np.trapz(y=result_f1[:,14], x=T1))
AUC2 = (np.trapz(y=result_f2[:,13], x=T1))/(np.trapz(y=result_f2[:,14], x=T1))
AUC3 = (np.trapz(y=result_f3[:,13], x=T1))/(np.trapz(y=result_f3[:,14], x=T1))
AUC4 = (np.trapz(y=result_f4[:,13], x=T1))/(np.trapz(y=result_f4[:,14], x=T1))
AUC5 = (np.trapz(y=result_f5[:,13], x=T1))/(np.trapz(y=result_f5[:,14], x=T1))
AUC6 = (np.trapz(y=result_f6[:,13], x=T1))/(np.trapz(y=result_f6[:,14], x=T1))
AUC7 = (np.trapz(y=result_f7[:,13], x=T1))/(np.trapz(y=result_f7[:,14], x=T1))
AUC8 = (np.trapz(y=result_f8[:,13], x=T1))/(np.trapz(y=result_f8[:,14], x=T1))
AUC9 = (np.trapz(y=result_f9[:,13], x=T1))/(np.trapz(y=result_f9[:,14], x=T1))

AUCratio = [AUC1,AUC2, AUC3,AUC4,AUC5,AUC6,AUC7,AUC8,AUC9]


plt.rcParams.update({'font.size': 25})
r1 = np.arange(9)
prop_iter = iter(plt.rcParams['axes.prop_cycle'])
for i in range(0,len(r1)):
  plt.bar(r1[i],AUCratio[i],color=next(prop_iter)['color'])

plt.xticks([r for r in range(9)], ['1','2','3','4', '5', '6', '7', '8', '9'])
plt.ylabel('Osteoblast/Osteoclast AUC')
plt.show()


# Pre-Osteoblast/Osteoblast area under curve
AUC1 = (np.trapz(y=result_f1[:,12], x=T1))/(np.trapz(y=result_f1[:,13], x=T1))
AUC2 = (np.trapz(y=result_f2[:,12], x=T1))/(np.trapz(y=result_f2[:,13], x=T1))
AUC3 = (np.trapz(y=result_f3[:,12], x=T1))/(np.trapz(y=result_f3[:,13], x=T1))
AUC4 = (np.trapz(y=result_f4[:,12], x=T1))/(np.trapz(y=result_f4[:,13], x=T1))
AUC5 = (np.trapz(y=result_f5[:,12], x=T1))/(np.trapz(y=result_f5[:,13], x=T1))
AUC6 = (np.trapz(y=result_f6[:,12], x=T1))/(np.trapz(y=result_f6[:,13], x=T1))
AUC7 = (np.trapz(y=result_f7[:,12], x=T1))/(np.trapz(y=result_f7[:,13], x=T1))
AUC8 = (np.trapz(y=result_f8[:,12], x=T1))/(np.trapz(y=result_f8[:,13], x=T1))
AUC9 = (np.trapz(y=result_f9[:,12], x=T1))/(np.trapz(y=result_f9[:,13], x=T1))

AUCratio = [AUC1,AUC2, AUC3,AUC4,AUC5,AUC6,AUC7,AUC8,AUC9]


plt.rcParams.update({'font.size': 25})
r1 = np.arange(9)
prop_iter = iter(plt.rcParams['axes.prop_cycle'])
for i in range(0,len(r1)):
  plt.bar(r1[i],AUCratio[i],color=next(prop_iter)['color'])

plt.xticks([r for r in range(9)], ['1','2','3','4', '5', '6', '7', '8', '9'])
plt.ylabel('Pre-Osteoblast/Osteoblast AUC')
plt.show()

# change in osteoblast from baseline
AUC1 = ((np.trapz(y=result_f1[:,13], x=T1))-(np.trapz(y=result_f1[:,13], x=T1)))/(np.trapz(y=result_f1[:,13], x=T1))
AUC2 = ((np.trapz(y=result_f2[:,13], x=T1))-(np.trapz(y=result_f1[:,13], x=T1)))/(np.trapz(y=result_f1[:,13], x=T1))
AUC3 = ((np.trapz(y=result_f3[:,13], x=T1))-(np.trapz(y=result_f1[:,13], x=T1)))/(np.trapz(y=result_f1[:,13], x=T1))
AUC4 = ((np.trapz(y=result_f4[:,13], x=T1))-(np.trapz(y=result_f1[:,13], x=T1)))/(np.trapz(y=result_f1[:,13], x=T1))
AUC5 = ((np.trapz(y=result_f5[:,13], x=T1))-(np.trapz(y=result_f1[:,13], x=T1)))/(np.trapz(y=result_f1[:,13], x=T1))
AUC6 = ((np.trapz(y=result_f6[:,13], x=T1))-(np.trapz(y=result_f1[:,13], x=T1)))/(np.trapz(y=result_f1[:,13], x=T1))
AUC7 = ((np.trapz(y=result_f7[:,13], x=T1))-(np.trapz(y=result_f1[:,13], x=T1)))/(np.trapz(y=result_f1[:,13], x=T1))
AUC8 = ((np.trapz(y=result_f8[:,13], x=T1))-(np.trapz(y=result_f1[:,13], x=T1)))/(np.trapz(y=result_f1[:,13], x=T1))
AUC9 = ((np.trapz(y=result_f9[:,13], x=T1))-(np.trapz(y=result_f1[:,13], x=T1)))/(np.trapz(y=result_f1[:,13], x=T1))

AUCratio = [AUC1,AUC2, AUC3,AUC4,AUC5,AUC6,AUC7,AUC8,AUC9]

plt.rcParams.update({'font.size': 25})
r1 = np.arange(9)
prop_iter = iter(plt.rcParams['axes.prop_cycle'])
for i in range(0,len(r1)):
  plt.bar(r1[i],AUCratio[i]*100,color=next(prop_iter)['color'])

plt.xticks([r for r in range(9)], ['1','2','3','4', '5', '6', '7', '8', '9'])
plt.ylabel('Change in Osteoblast AUC')
plt.show()


# change in osteoclast from baseline
AUC1 = ((np.trapz(y=result_f1[:,14], x=T1))-(np.trapz(y=result_f1[:,14], x=T1)))/(np.trapz(y=result_f1[:,14], x=T1))
AUC2 = ((np.trapz(y=result_f2[:,14], x=T1))-(np.trapz(y=result_f1[:,14], x=T1)))/(np.trapz(y=result_f1[:,14], x=T1))
AUC3 = ((np.trapz(y=result_f3[:,14], x=T1))-(np.trapz(y=result_f1[:,14], x=T1)))/(np.trapz(y=result_f1[:,14], x=T1))
AUC4 = ((np.trapz(y=result_f4[:,14], x=T1))-(np.trapz(y=result_f1[:,14], x=T1)))/(np.trapz(y=result_f1[:,14], x=T1))
AUC5 = ((np.trapz(y=result_f5[:,14], x=T1))-(np.trapz(y=result_f1[:,14], x=T1)))/(np.trapz(y=result_f1[:,14], x=T1))
AUC6 = ((np.trapz(y=result_f6[:,14], x=T1))-(np.trapz(y=result_f1[:,14], x=T1)))/(np.trapz(y=result_f1[:,14], x=T1))
AUC7 = ((np.trapz(y=result_f7[:,14], x=T1))-(np.trapz(y=result_f1[:,14], x=T1)))/(np.trapz(y=result_f1[:,14], x=T1))
AUC8 = ((np.trapz(y=result_f8[:,14], x=T1))-(np.trapz(y=result_f1[:,14], x=T1)))/(np.trapz(y=result_f1[:,14], x=T1))
AUC9 = ((np.trapz(y=result_f9[:,14], x=T1))-(np.trapz(y=result_f1[:,14], x=T1)))/(np.trapz(y=result_f1[:,14], x=T1))

AUCratio = [AUC1,AUC2, AUC3,AUC4,AUC5,AUC6,AUC7,AUC8,AUC9]


plt.rcParams.update({'font.size': 25})
r1 = np.arange(9)
prop_iter = iter(plt.rcParams['axes.prop_cycle'])
for i in range(0,len(r1)):
  plt.bar(r1[i],AUCratio[i]*100,color=next(prop_iter)['color'])

plt.xticks([r for r in range(9)], ['1','2','3','4', '5', '6', '7', '8', '9'])
#plt.ylabel('Osteoblast/Osteoclast AUC')
#plt.ylabel('Pre-Osteoblast/Osteoblast AUC')
#plt.ylabel('Change in Osteoblast AUC')
plt.ylabel('Change in Osteoclast AUC')
plt.show()


#input('')

plt.rcParams.update({'font.size': 25})
#plt.xticks(np.arange(0, 30, step=7))

BFRC = []
BFRC1 = []
BFRC2 = []

for i in range(len(T1)-1):
    print(i)
    BFR = (result_f4[i+1,15]-result_f4[i,15])/(result_f1[i+1,15]-result_f1[i,15])
    BFRC.append(BFR)
    BFRC1.append((result_f4[i+1,15]-result_f4[i,15])/0.001)
    BFRC2.append((result_f1[i+1,15]-result_f1[i,15])/0.001)


#plt.plot(T1[0:-1], BFRC, linewidth=3)
plt.plot(T1[0:-1], BFRC1,T1[0:-1], BFRC2, linewidth=3)
plt.xlabel('TIME (days)')
plt.ylabel('Change in BFR')
plt.legend(['without butyrate', 'With butyrate + B'])
plt.show()

plt.xticks(np.arange(0, 30, step=7))
plt.plot(T1, result_f1[:,15],T1, result_f2[:,15],T1, result_f3[:,15], T1, result_f4[:,15],T1, result_f5[:,15],T1, result_f6[:,15],T1, result_f7[:,15],T1, result_f8[:,15],T1, result_f9[:,15], linewidth=3)
#plt.plot(T1, result_f1[:,15],T1, result_f2[:,15],T1, result_f3[:,15], T1, result_f4[:,15],T1, result_f5[:,15],T1, result_f6[:,15],T1, result_f7[:,15],T1, result_f8[:,15],T1, result_f9[:,15], linewidth=20)
plt.xlabel('Time (days)')
plt.ylabel('Relative bone volume (%)')
#plt.grid()
#plt.legend(['without butyrate','with butyrate', 'osteoblast proliferation', 'osteoblast differentiation', 'osteoblast differentiation (wnt)', 'osteoclast death', 'osteoclast differentiation','bone formation','bone resorption'],loc='upper left')
#plt.legend(['Case 1: without butyrate','Case 2: with butyrate', 'Case 3: osteoblast proliferation', 'Case 4: osteoblast differentiation', 'Case 5: osteoblast differentiation via Wnto10b', 'Case 6: osteoclast death', 'Case 7: osteoclast differentiation','Case 8: bone formation','Case 9: bone resorption'],loc='upper left')

plt.errorbar([28],[136.6], yerr=[6.47], fmt='o',color='r', elinewidth=3, markersize=10, capsize=6, capthick=3, barsabove= False)
#plt.ylim([80,400])
plt.show()

plt.xticks(np.arange(0, 30, step=7))
plt.plot(T1, result_f1[:,11],T1, result_f2[:,11],T1, result_f3[:,11], T1, result_f4[:,11],T1, result_f5[:,11],T1, result_f6[:,11],T1, result_f7[:,11],T1, result_f8[:,11],T1, result_f9[:,11], linewidth=3)
plt.xlabel('Time (days)')
plt.ylabel('Osteocyte cells')
#plt.legend(['without butyrate','with butyrate', 'osteoblast proliferation', 'osteoblast differentiation', 'osteoblast differentiation (wnt)', 'osteoclast death', 'osteoclast differentiation','bone formation','bone resorption'],loc='center right')
plt.show()

plt.xticks(np.arange(0, 30, step=7))
plt.plot(T1, result_f1[:,12],T1, result_f2[:,12],T1, result_f3[:,12], T1, result_f4[:,12],T1, result_f5[:,12],T1, result_f6[:,12],T1, result_f7[:,12],T1, result_f8[:,12],T1, result_f9[:,12], linewidth=3)
plt.xlabel('Time (days)')
plt.ylabel('Pre-osteoblast cells')
#plt.legend(['without butyrate','with butyrate', 'osteoblast proliferation', 'osteoblast differentiation', 'osteoblast differentiation (wnt)', 'osteoclast death', 'osteoclast differentiation','bone formation','bone resorption'],loc='center right')
plt.show()

plt.xticks(np.arange(0, 30, step=7))
plt.plot(T1, result_f1[:,13],T1, result_f2[:,13],T1, result_f3[:,13], T1, result_f4[:,13],T1, result_f5[:,13],T1, result_f6[:,13],T1, result_f7[:,13],T1, result_f8[:,13],T1, result_f9[:,13], linewidth=3)
plt.xlabel('Time (days)')
plt.ylabel('Osteoblast cells')
#plt.legend(['without butyrate','with butyrate', 'osteoblast proliferation', 'osteoblast differentiation', 'osteoblast differentiation (wnt)', 'osteoclast death', 'osteoclast differentiation','bone formation','bone resorption'],loc='center right')
#plt.legend(['Without butyrate','With butyrate', 'With butyrate + O', 'With butyrate + B'], loc='center right')
plt.show()

plt.xticks(np.arange(0, 30, step=7))
plt.plot(T1, result_f1[:,14],T1, result_f2[:,14],T1, result_f3[:,14], T1, result_f4[:,14],T1, result_f5[:,14],T1, result_f6[:,14],T1, result_f7[:,14],T1, result_f8[:,14],T1, result_f9[:,14], linewidth=3)
plt.xlabel('Time (days)')
plt.ylabel('Osteoclast cells')
#plt.legend(['without butyrate','with butyrate', 'osteoblast proliferation', 'osteoblast differentiation', 'osteoblast differentiation (wnt)', 'osteoclast death', 'osteoclast differentiation','bone formation','bone resorption'],loc='center right')
#plt.legend(['Without butyrate','With butyrate', 'With butyrate + O', 'With butyrate + B'], loc='center right')
plt.show()

plt.xticks(np.arange(0, 30, step=7))
plt.plot(T1, result_f1[:,10],'--',color = 'magenta', linewidth=3)
plt.plot(T1, result_f2[:,10],color = 'purple',  linewidth=3)
plt.xlabel('Time (days)')
plt.ylabel('Wnt10b fold change')
plt.legend(['without butyrate','with butyrate'])
plt.show()