import numpy as np
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import random,pickle

result1 = np.array(pickle.load(open('result1a.p', 'rb')))
result2 = np.array(pickle.load(open('result2a.p', 'rb')))
result3 = np.array(pickle.load(open('result3a.p', 'rb')))
result4 = np.array(pickle.load(open('result4a.p', 'rb')))
result5 = np.array(pickle.load(open('result5a.p', 'rb')))
result6 = np.array(pickle.load(open('result6a.p', 'rb')))
T = np.array(pickle.load(open('time1a.p', 'rb')))

plt.rcParams.update({'font.size': 25,'legend.fontsize': 15})
plt.yticks(np.arange(42, 60, step=2.5))
plt.plot(T, (result1[:,8]/(result1[:,5] + result1[:,8]))*100,T, (result2[:,8]/(result2[:,5] + result2[:,8]))*100,T, (result3[:,8]/(result3[:,5] + result3[:,8]))*100,T, (result4[:,8]/(result4[:,5] + result4[:,8]))*100,T, (result5[:,8]/(result5[:,5] + result5[:,8]))*100,T, (result6[:,8]/(result6[:,5] + result6[:,8]))*100, linewidth=3)
plt.legend(['0.18 \u03bcM','0.25 \u03bcM','0.30 \u03bcM','0.35 \u03bcM','0.40 \u03bcM','0.50 \u03bcM'])
plt.xlabel('Time (days)')
plt.ylabel('Bone Tregs (%)')
plt.savefig("images/IECR/" + "12.png", dpi = 300, bbox_inches='tight')
plt.show()

result1 = np.array(pickle.load(open('result7a.p', 'rb')))
result2 = np.array(pickle.load(open('result8a.p', 'rb')))
result3 = np.array(pickle.load(open('result9a.p', 'rb')))
result4 = np.array(pickle.load(open('result10a.p', 'rb')))
result5 = np.array(pickle.load(open('result11a.p', 'rb')))
result6 = np.array(pickle.load(open('result12a.p', 'rb')))

plt.rcParams.update({'font.size': 25,'legend.fontsize': 15})
plt.yticks(np.arange(22, 44, step=2))
plt.plot(T, (result1[:,8]/(result1[:,5] + result1[:,8]))*100,T, (result2[:,8]/(result2[:,5] + result2[:,8]))*100,T, (result3[:,8]/(result3[:,5] + result3[:,8]))*100,T, (result4[:,8]/(result4[:,5] + result4[:,8]))*100,T, (result5[:,8]/(result5[:,5] + result5[:,8]))*100,T, (result6[:,8]/(result6[:,5] + result6[:,8]))*100, linewidth=3)
plt.legend(['-0.11 \u03bcM','-0.13 \u03bcM','-0.15 \u03bcM','-0.17 \u03bcM','-0.19 \u03bcM','-0.21 \u03bcM'])
plt.xlabel('Time (days)')
plt.ylabel('Bone Tregs (%)')
plt.savefig("images/IECR/" + "13.png", dpi = 300, bbox_inches='tight')
plt.show()