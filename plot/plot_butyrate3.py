import numpy as np
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import random,pickle

result1 = np.array(pickle.load(open('resultP.p', 'rb')))
result2 = np.array(pickle.load(open('resultN.p', 'rb')))
T = np.array(pickle.load(open('time1.p', 'rb')))

plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
plt.xticks(np.arange(0, 30, step=7))
plt.yticks(np.arange(0, 3.6, step=0.5))
plt.plot(T, result1[:,9],T, result2[:,9], linewidth=3)
plt.legend(['0.18 \u03bcM','-0.11 \u03bcM'])
plt.ylabel('TGF-β fold change')
plt.xlabel('Time (days)')
plt.ylim([0,3.5])
plt.savefig("images/IECR/" + "14.png", dpi = 300, bbox_inches='tight')
plt.show()

plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
plt.xticks(np.arange(0, 30, step=7))
plt.yticks(np.arange(0, 3.6, step=0.5))
plt.plot(T, result1[:,10],T, result2[:,10], linewidth=3)
plt.legend(['0.18 \u03bcM','-0.11 \u03bcM'])
plt.ylabel('Wnt10b fold change')
plt.xlabel('Time (days)')
plt.ylim([0,3.5])
plt.savefig("images/IECR/" + "15.png", dpi = 300, bbox_inches='tight')
plt.show()