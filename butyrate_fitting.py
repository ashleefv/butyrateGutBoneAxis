import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import odeint
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


tspan = [1, 2, 3]
C_data = [1.46, 21.50, 27.45]      # data without butyrate
#C_data = [1.46, 26.35, 46.40]      # data with butyrate

# without butyrate = 0.1717933, with butyrate = 0.30026203, combined = 0.12846854
def fitfunc(t, k):
    def myode(x, t):
        y = [0.0 for i in range(len(x))]
        y[0] = k * (100 - x[0])
        # y[0] = 0.1717933 * (100 - x[0]) + k * (100 - x[0])
        return y

    x = (C_data[0])     #initial value
    C_sol = odeint(myode, x, t)
    return C_sol[:,0]

k_fit, kcov = curve_fit(fitfunc, tspan, C_data, p0=0)
print (k_fit, kcov)

tfit = np.linspace(1,3,100)
fit = fitfunc(tfit, k_fit[0])

plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
plt.plot(tspan, C_data, 'ro', label='data', markersize=12)
plt.plot(tfit, fit, 'b-', label='fit', linewidth=3)
plt.legend(loc='upper left')
plt.xlabel('Time (days)')
plt.ylabel('Tregs (%)')
plt.ylim([0,50])
#plt.savefig("images/IECR/" + "17.png", dpi = 300, bbox_inches='tight')
plt.show()