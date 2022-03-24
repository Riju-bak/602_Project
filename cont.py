# Kinetic Parameters
KiC = 1.15
KiN = 0.073
KiX = 152
KiS = 100
KN = 0.033
KO1 = 0.65
KO2 = 5.50
KS = 0.077
KSE = 0.0032
KSL = 0.021
mS = 0.012
rE = 0.31
rL = 0.78
YCS = 0.89
YLS = 0.47
YXN = 27.0
YXS = 1.52
aE = 0.002
aL = 0.017
bCmax = 0.090
mmax = 0.26

import numpy as np
from scipy import stats
from scipy import optimize
from pickle import TRUE
import pandas as pd

# Operating Conditions Input
n = int(input("Number of stage: "))
O = 25

# n-stage
data = []
for a in range (0, n+1):
    if a == 0:
        SF = 0
        D = 0
        V = 0
        S = 0
        N = 0
        Xf = 0
        L = 0
        E = 0
        C = 0 
        F = 0
    else:
        print("Inputs for stage {}".format(a))
        SF = float(input("Glucose Concentration in Feed {} (g/L): ".format(a)))
        D = float(input("Dilution Rate of Stage {}: ".format(a)))
        if a == 1:
            V = float(input("Volume of stage-1 reactor (L): "))
            CN = float(input("Ratio C/N in Feed 1: "))
            NF = (SF*6*12/180) / CN
        else:
            VRatio = float(input("Volume Ratio V{}/V{}: ".format(a, a-1)))
            V = V*VRatio
            F = V*D
    name = ['SF', 'D', 'V', 'S', 'N', 'Xf', 'L', 'E', 'C']
    value = [SF, D, V, S, N, Xf, L, E, C]
    data.append((dict(zip(name, value))))
df = pd.DataFrame(data)

     #Solve Non-linear Equations
for a in range (1, n+1, 1):     
    def equations(vars):
        FB, S, N, Xf, L, E, C  = vars
        m = mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*(KiS/(KiS+S))*(KiX/(KiX+Xf))
        bLC = bCmax*(KiN/(KiN+N))*(O/(KO2+O))*(S/(KS+S))*(KiS/(KiS+S))*(KiX/(KiX+Xf))*((KiC - (C/Xf))/KiC)
        bC = (1-rL)*bLC
        bL = rL*bLC - KSL*(L/(Xf+L))*(O/(KO2+O))
        bE = rE*bLC*(O/(KO2+O)) - KSE*(E/(Xf+L))*(O/(KO2+O))
        qS = m/YXS + (O/(KO1*Xf + O))*(S/(KS + S))*mS + bC/YCS + (aL*m + bL)/YLS
        qN = m/YXN
        qL = aL*m + bL
        qE = aE*m + bE
        qC = 1.88*(1-rL)*bLC
        FS = df.loc[(a), 'D']*df.loc[(a), 'V'] - FB
        F = df.loc[(a), 'D']*df.loc[(a), 'V']
        eq1 = FB - ((df.loc[(a), 'V']/1000) * ((7.14/YXN)*m*Xf + 1.59*bC*Xf))
        eq2 = F*df.loc[(a-1), 'Xf'] - F*Xf + m*Xf*df.loc[(a), 'V']
        eq3 = F*df.loc[(a-1), 'S'] + FS*df.loc[(a), 'SF'] - F*S - qS*Xf*df.loc[(a), 'V']
        eq4 = F*df.loc[(a-1), 'N'] + FS*NF - F*N - qN*Xf*df.loc[(a), 'V']
        eq5 = F*df.loc[(a-1), 'L'] - F*L + qL*Xf*df.loc[(a), 'V']
        eq6 = F*df.loc[(a-1), 'E'] - F*E + qE*Xf*df.loc[(a), 'V']
        eq7 = F*df.loc[(a-1), 'C'] - F*C + qC*Xf*df.loc[(a), 'V']
        return (eq1, eq2, eq3, eq4, eq5, eq6, eq7)
    def f(vars):
        return abs(np.array(sum(equations(vars))**2)-0)
    FB, S, N, Xf, L, E, C =  optimize.fsolve(equations, (optimize.fmin(f, (0.02, 0.095, 0.01, 150, 100, 56, 100), xtol=0.01, maxiter=1000)))
    df.loc[a, 'S'], df.loc[a, 'N'], df.loc[a, 'Xf'], df.loc[a, 'L'], df.loc[a, 'E'], df.loc[a, 'C'] = S, N, Xf, L, E, C
df = np.round(df, decimals = 2)
print(df)

#Evaluation Results 
t = (df['V'].sum()) / (D*V)
EPA_titer = E*100/(Xf+L)
EPA_rate = E/t

print("EPA titer: ", round(EPA_titer))
print("EPA_rate: ", round(EPA_rate))
