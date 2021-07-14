import numpy as np
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import random,pickle


# for positive butyrate dose
result = np.array(pickle.load(open('result.p', 'rb')))               # for input bI1 = 0.18
result1 = np.array(pickle.load(open('result1.p', 'rb')))             # for input bI1 = 0.25
result2 = np.array(pickle.load(open('result2.p', 'rb')))             # for input bI1 = 0.11
result3 = np.array(pickle.load(open('resulth.p', 'rb')))             # for input bI1 = 0
T = np.array(pickle.load(open('time1.p', 'rb')))


# butyrate in intestine
plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
plt.xticks(np.arange(0, 30, step=7))
plt.plot(T, result[:,0], color = 'orange',linewidth=3)
plt.fill_between(T, result1[:,0], result2[:,0], color = 'k', alpha = 0.1)
plt.legend(['butyrate in intestine'], loc='lower right')
plt.xlabel('Time (days)')
plt.ylabel('$\Delta$ Butyrate \u03bcM')
plt.ylim([0,0.42])
plt.errorbar([28], [0.18], yerr=[0.07], fmt='s',color='r', fillstyle = 'none', elinewidth=3, markersize=10, capsize=6, capthick=3, barsabove= False)
#plt.savefig("images/IECR/" + "1.png", dpi = 300, bbox_inches='tight')
plt.show()

# butyrate in blood
plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
plt.xticks(np.arange(0, 30, step=7))
plt.plot(T, result[:,1], color = 'orange',linewidth=3)
plt.fill_between(T, result1[:,1], result2[:,1], color = 'k', alpha = 0.1)
plt.legend(['butyrate in blood'], loc='lower right')
plt.xlabel('Time (days)')
plt.ylabel('$\Delta$ Butyrate \u03bcM')
plt.ylim([0,0.42])
plt.errorbar([28], [0.29], yerr=[0.09], fmt='s',color='r', fillstyle = 'none', elinewidth=3, markersize=10, capsize=6, capthick=3, barsabove= False)
#plt.savefig("images/IECR/" + "2.png", dpi = 300, bbox_inches='tight')
plt.show()

# butyrate in bone
plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
plt.xticks(np.arange(0, 30, step=7))
plt.plot(T, result[:,2], color = 'orange',linewidth=3)
plt.fill_between(T, result1[:,2], result2[:,2], color = 'k', alpha = 0.1)
plt.legend(['butyrate in bone'], loc='lower right')
plt.xlabel('Time (days)')
plt.ylabel('$\Delta$ Butyrate \u03bcM')
plt.ylim([0,0.42])
plt.errorbar([28], [0.29], yerr=[0], fmt='s',color='r', fillstyle = 'none', elinewidth=3, markersize=10, capsize=0, capthick=3, barsabove= False)
#plt.savefig("images/IECR/" + "3.png", dpi = 300, bbox_inches='tight')
plt.show()

# Tregs in blood
plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
plt.xticks(np.arange(0, 30, step=7))
plt.yticks(np.arange(13, 18, step=1))
plt.plot(T, (result3[:,7]/(result3[:,4] + result3[:,7]))*100 ,T, (result[:,7]/(result[:,4] + result[:,7]))*100,linewidth=3)
plt.fill_between(T, (result1[:,7]/(result1[:,4] + result1[:,7]))*100, (result2[:,7]/(result2[:,4] + result2[:,7]))*100, color = 'k', alpha = 0.1)
plt.errorbar([28, 28],[13, 16], yerr=[0.6,0.4], fmt='s', fillstyle = 'none', color='r', elinewidth=3, markersize=10, capsize=6, capthick=3, barsabove= False)
plt.legend(['Control', 'Control + butyrate'], loc='upper left')
plt.xlabel('Time (days)')
plt.ylabel('Blood Tregs (%)')
#plt.savefig("images/IECR/" + "4.png", dpi = 300, bbox_inches='tight')
plt.show()

# Tregs in bone
plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
plt.xticks(np.arange(0, 30, step=7))
plt.plot(T, (result3[:,8]/(result3[:,5] + result3[:,8]))*100 ,T, (result[:,8]/(result[:,5] + result[:,8]))*100,linewidth=3)
plt.fill_between(T, (result1[:,8]/(result1[:,5] + result1[:,8]))*100, (result2[:,8]/(result2[:,5] + result2[:,8]))*100, color = 'k', alpha = 0.1)
plt.errorbar([28],[42], yerr=[1], fmt='s',color='r', fillstyle = 'none', elinewidth=3, markersize=10, capsize=6, capthick=3, barsabove= False)
plt.errorbar([28],[48], yerr=[2], fmt='o',color='r', elinewidth=3, markersize=10, capsize=6, capthick=3, barsabove= False)
plt.legend(['Control', 'Control + butyrate'], loc='upper left')
plt.xlabel('Time (days)')
plt.ylabel('Bone Tregs (%)')
#plt.savefig("images/IECR/" + "5.png", dpi = 300, bbox_inches='tight')
plt.show()


# for negative butyrate dose

result0 = np.array(pickle.load(open('result0.p', 'rb')))             # for input bI1 = -0.11
result01 = np.array(pickle.load(open('result01.p', 'rb')))           # for input bI1 = -0.13
result02 = np.array(pickle.load(open('result02.p', 'rb')))           # for input bI1 = -0.09
result3 = np.array(pickle.load(open('resulth.p', 'rb')))             # for input bI1 = 0
T = np.array(pickle.load(open('time1.p', 'rb')))

# butyrate in intestine
plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
plt.xticks(np.arange(0, 30, step=7))
plt.plot(T, result0[:,0], color = 'orange',linewidth=3)
plt.fill_between(T, result01[:,0], result02[:,0], color = 'k', alpha = 0.1)
plt.legend(['butyrate in intestine'], loc='upper right')
plt.xlabel('Time (days)')
plt.ylabel('$\Delta$ Butyrate \u03bcM')
plt.ylim([-0.22,0.01])
#plt.savefig("images/IECR/" + "6.png", dpi = 300, bbox_inches='tight')
plt.show()

# butyrate in blood
plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
plt.xticks(np.arange(0, 30, step=7))
plt.plot(T, result0[:,1], color = 'orange',linewidth=3)
plt.fill_between(T, result01[:,1], result02[:,1], color = 'k', alpha = 0.1)
plt.legend(['butyrate in blood'], loc='upper right')
plt.xlabel('Time (days)')
plt.ylabel('$\Delta$ Butyrate \u03bcM')
plt.ylim([-0.22,0.01])
plt.errorbar([28], [-0.17722], yerr=[0.032], fmt='s',color='r', fillstyle = 'none', elinewidth=3, markersize=10, capsize=6, capthick=3, barsabove= False)
#plt.savefig("images/IECR/" + "7.png", dpi = 300, bbox_inches='tight')
plt.show()

# butyrate in bone
plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
plt.xticks(np.arange(0, 30, step=7))
plt.plot(T, result0[:,2], color = 'orange',linewidth=3)
plt.fill_between(T, result01[:,2], result02[:,2], color = 'k', alpha = 0.1)
plt.legend(['butyrate in bone'], loc='upper right')
plt.xlabel('Time (days)')
plt.ylabel('$\Delta$ Butyrate \u03bcM')
plt.ylim([-0.22,0.01])
#plt.savefig("images/IECR/" + "8.png", dpi = 300, bbox_inches='tight')
plt.show()


# Tregs in bone
plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
plt.xticks(np.arange(0, 30, step=7))
plt.yticks(np.arange(37, 44, step=1))
plt.plot(T, (result3[:,8]/(result3[:,5] + result3[:,8]))*100 ,T, (result0[:,8]/(result0[:,5] + result0[:,8]))*100,linewidth=3)
plt.fill_between(T, (result01[:,8]/(result01[:,5] + result01[:,8]))*100, (result02[:,8]/(result02[:,5] + result02[:,8]))*100, color = 'k', alpha = 0.1)
plt.errorbar([28],[42], yerr=[1], fmt='s',color='r', fillstyle = 'none', elinewidth=3, markersize=10, capsize=6, capthick=3, barsabove= False)
plt.errorbar([28],[37.5], yerr=[1], fmt='o',color='r', elinewidth=3, markersize=10, capsize=6, capthick=3, barsabove= False)
plt.legend(['Control', 'Control + butyrate'], loc='lower left')
plt.xlabel('Time (days)')
plt.ylabel('Bone Tregs (%)')
#plt.savefig("images/IECR/" + "9.png", dpi = 300, bbox_inches='tight')
plt.show()