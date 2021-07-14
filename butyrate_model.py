import numpy as np
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pickle


from butyrate_model_constants import *


result = []
output = []
for i in range(Nc):

    def diff(x,T):

        # Differential equations

        y = [0.0 for i in range(len(x))]

        # butyrate in intestine
        y[0] = FB1 - muB * x[0] - AB12 * x[0]

        # butyrate in blood
        y[1] = AB12 * x[0] - muB * x[1] - AB23 * x[1]

        # butyrate in bone
        y[2] = AB23 * x[1] - muB * x[2]

        # naive T cells in intestine
        y[3] = FN1_minus - b_minus * gamma * x[3] - muN * x[3] - b_plus * x[3] * x[0] + FN1_plus

        # naive T cells in Blood
        y[4] = FN2_minus - b_minus * gamma * x[4] - muN * x[4] - b_plus * x[4] * x[1] + FN2_plus

        # naive T cells in Bone
        y[5] = FN3_minus - b_minus * gamma * x[5] - muN * x[5] - b_plus * x[5] * x[2] + FN3_plus

        # regulatory T cells in intestine
        y[6] = b_minus * gamma * x[3] + b_plus * x[3] * x[0] - deltaT12 * x[6] - muT * x[6]

        # regulatory T cells in blood
        y[7] = deltaT12 * x[6] - muT * x[7] - deltaT23 * x[7] + b_minus * gamma * x[4] + b_plus * x[4] * x[1]

        # regulatory T cells in bone
        y[8] = deltaT23 * x[7] - muT * x[8] + b_plus * x[5] * x[2] + b_minus * gamma * x[5]

        # TGF-beta production and decay
        y[9] = muT_beta * (1 - x[9]) + VT_beta * (x[8] / x8 - 1) / (kT_beta - (x[8] / x8 - 1))

        # Wnt10b formation
        y[10] = muW * (1 - x[10]) + rhoW * (x[9] - 1)


        # S Osteocytes
        y[11] = (alpha_1*math.pow(x[13],g_31)*max((1-x[11]/K_S),0))*100/sc

        # P Pre-osteoblast
        # butyrate increase pre-osteoblast to osteoblast differentiation
        if case == 4:
            y[12] = (alpha_2*math.pow(x[11],g_21)*math.pow(max((1-x[11]/K_S),0),g_22) + alpha_3*math.pow(x[12],g_32)*max((1-x[11]/K_S),0) + (kg_2*((x[10]-x10)/((x[10]-x10) + kg_4)))*math.pow(x[12],f_12)  - (beta_1 + x[2] * kod)*math.pow(x[12],f_12)*math.pow(x[14],f_14) - (kg_1*((x[10]-x10)/((x[10]-x10) + kg_4)))*math.pow(x[12],f_12) - delta*x[12])*100/sc

        # butyrate increase wnt10b dependent pre-osteoblast to osteoblast differentiation
        elif case == 5:
            y[12] = (alpha_2*math.pow(x[11],g_21)*math.pow(max((1-x[11]/K_S),0),g_22) + alpha_3*math.pow(x[12],g_32)*max((1-x[11]/K_S),0) + (kg_2*((x[10]-x10)/((x[10]-x10) + kg_4)))*math.pow(x[12],f_12)  - beta_1*math.pow(x[12],f_12)*math.pow(x[14],f_14) - ((kg_1*((x[10]-x10)/((x[10]-x10) + kg_4))) + x[2] * kodw)*math.pow(x[12],f_12) - delta*x[12])*100/sc

        # default equation
        else:
            y[12] = (alpha_2*math.pow(x[11],g_21)*math.pow(max((1-x[11]/K_S),0),g_22) + alpha_3*math.pow(x[12],g_32)*max((1-x[11]/K_S),0) + (kg_2*((x[10]-x10)/((x[10]-x10) + kg_4)))*math.pow(x[12],f_12)  - beta_1*math.pow(x[12],f_12)*math.pow(x[14],f_14) - (kg_1*((x[10]-x10)/((x[10]-x10) + kg_4)))*math.pow(x[12],f_12) - delta*x[12])*100/sc

        # O Osteoblast
        # butyrate proliferate osteoblast
        if case == 3:
            y[13] = (kop*x[2]*x[13] - alpha_1*math.pow(x[13],g_31)*max((1-x[11]/K_S),0) + beta_1*math.pow(x[12],f_12)*math.pow(x[14],f_14) + (kg_1*((x[10]-x10)/((x[10]-x10) + kg_4)))*math.pow(x[12],f_12) - (beta_2-(kg_3*((x[10]-x10)/((x[10]-x10) + kg_4))))*math.pow(x[13],f_23))*100/sc

        # butyrate increase pre-osteoblast to osteoblast differentiation
        elif case == 4:
            y[13] = (- alpha_1*math.pow(x[13],g_31)*max((1-x[11]/K_S),0) + (beta_1 + x[2] * kod)*math.pow(x[12],f_12)*math.pow(x[14],f_14) + (kg_1*((x[10]-x10)/((x[10]-x10) + kg_4)))*math.pow(x[12],f_12) - (beta_2-(kg_3*((x[10]-x10)/((x[10]-x10) + kg_4))))*math.pow(x[13],f_23))*100/sc

        # butyrate increase wnt10b dependent pre-osteoblast to osteoblast differentiation
        elif case == 5:
            y[13] = (- alpha_1*math.pow(x[13],g_31)*max((1-x[11]/K_S),0) + beta_1*math.pow(x[12],f_12)*math.pow(x[14],f_14) + ((kg_1*((x[10]-x10)/((x[10]-x10) + kg_4))) + x[2] * kodw)*math.pow(x[12],f_12) - (beta_2-(kg_3*((x[10]-x10)/((x[10]-x10) + kg_4))))*math.pow(x[13],f_23))*100/sc

        # default equation
        else:
            y[13] = (- alpha_1*math.pow(x[13],g_31)*max((1-x[11]/K_S),0) + beta_1*math.pow(x[12],f_12)*math.pow(x[14],f_14) + (kg_1*((x[10]-x10)/((x[10]-x10) + kg_4)))*math.pow(x[12],f_12) - (beta_2-(kg_3*((x[10]-x10)/((x[10]-x10) + kg_4))))*math.pow(x[13],f_23))*100/sc


        # C Osteoclast
        # butyrate inhibiting osteoclast by osteoclast death
        if case == 6:
            y[14] = (-kca*x[2]*x[14]+alpha_4*math.pow(x[11],g_41)*math.pow(x[12],g_42)*math.pow(epsilon+x[13],g_43)*math.pow(max((1-x[11]/K_S),0),g_44) - beta_3*math.pow(x[14],f_34))*100/sc

        # butyrate inhibiting osteoclast differentiation
        elif case == 7:
            y[14] = ((alpha_4-kcd*x[14])*math.pow(x[11],g_41)*math.pow(x[12],g_42)*math.pow(epsilon+x[13],g_43)*math.pow(max((1-x[11]/K_S),0),g_44) - beta_3*math.pow(x[14],f_34))*100/sc

        # default equation
        else:
            y[14] = (alpha_4*math.pow(x[11],g_41)*math.pow(x[12],g_42)*math.pow(epsilon+x[13],g_43)*math.pow(max((1-x[11]/K_S),0),g_44) - beta_3*math.pow(x[14],f_34))*100/sc


        # z bone volume
        # butyrate increase osteoblast mediated bone formation
        if case == 8:
            y[15] = (-k1 * x[14] + (k2 + x[2] * kzo) * x[13]) * 100 / sc

        # butyrate reduce osteoclast mediated bone resorption
        elif case == 9:
            y[15] = (-(k1 - x[2] * kzc) * x[14] + k2 * x[13]) * 100 / sc

        # default equation
        else:
            y[15] = (-k1 * x[14] + k2 * x[13]) * 100 / sc

        return y

    N = cyclelength       # Total number of days

    if i==0:
        x = (x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,S0,P0,B0,C0,z0)
    else:
        x = (result[-1,0], result[-1,1],result[-1,2],result[-1,3], result[-1,4],result[-1,5], result[-1,6],result[-1,7], result[-1,8],result[-1,9], result[-1,10], S0, result[-1,12] if result[-1,12]>=1 else 0,result[-1,13] if result[-1,13]>=1 else 0, result[-1,14] if result[-1,14]>=1 else 0,result[-1,15] if result[-1,15]>=1 else 0)


    T = np.arange(0.0, N, 0.001)

    result = np.array(odeint(diff, x, T))
    output.append(result)



# 1 cycle
if Nc==1:
    result_f = result

# 2 cycle
if Nc==2:
    result_f = np.concatenate((np.array(output[0]),np.array(output[1])))

# 6 cycle
if Nc==6:
    result_f = np.concatenate((np.array(output[0]),np.array(output[1]),np.array(output[2]),np.array(output[3]),np.array(output[4]),np.array(output[5])))

# 12 cycle
if Nc==12:
    result_f = np.concatenate((np.array(output[0]),np.array(output[1]),np.array(output[2]),np.array(output[3]),np.array(output[4]),np.array(output[5]), np.array(output[6]),np.array(output[7]),np.array(output[8]),np.array(output[9]),np.array(output[10]),np.array(output[11])))


#pickle.dump(result_f, open('case1.p', 'wb'))
T1 = np.arange(0.0, Nc*cyclelength, 0.001)

# plot
plt.rcParams.update({'font.size': 25})

plt.plot(T1, (result_f[:,8]/(result_f[:,5] + result_f[:,8]))*100, linewidth=3)
plt.xlabel('Time (days)')
plt.ylabel('Bone Tregs (%)')
plt.show()

plt.plot(T1, result_f[:,10], linewidth=3)
plt.xlabel('TIME (days)')
plt.ylabel('Wnt10b fold change')
plt.show()

plt.plot(T1, result_f[:,15], linewidth=3)
plt.xlabel('Time (days)')
plt.ylabel('Relative bone volume (%)')
plt.show()


