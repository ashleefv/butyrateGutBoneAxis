import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import odeint
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# given data we want to fit
#tspan = [1, 1.33, 1.48, 1.63, 1.81, 2, 2.2, 2.58, 2.77, 3]
#C_data = [1.46, 6.5, 9.8, 13.57, 16.6, 21.50, 22.93, 24.54, 26.34, 27.45]

#tspan = [1, 2, 3]
#C_data = [1.46, 21.50, 27.45]

tspan = [1, 2, 3]
C_data = [1.46, 26.35, 46.40]

# without butyrate = 0.1717933, with butyrate = 0.30026203, combined = 0.12846854
def fitfunc(t, k):
    def myode(x, t):
        y = [0.0 for i in range(len(x))]
        y[0] = k * (100-x[0])
        #y[0] = 0.1717933*(100-x[0]) + k * (100-x[0])
        return y

    x = (C_data[0])     #initial value
    C_sol = odeint(myode, x, t)
    return C_sol[:,0]

k_fit, kcov = curve_fit(fitfunc, tspan, C_data, p0=0)
print (k_fit, kcov)

tfit = np.linspace(1,3)
fit = fitfunc(tfit, k_fit[0])


plt.plot(tspan, C_data, 'ro', label='data')
plt.plot(tfit, fit, 'b-', label='fit')
plt.legend(loc='best')
plt.show()
