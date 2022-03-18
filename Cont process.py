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

# Number of Stage
#n = input("Number of stage: ")
#n = int(n)

# Operating Conditions Input
V = 1.5
O = 25
SF = 700
CN = 68
D = 0.009

#SF = input("Glucose Concentration in Feed 1: ")
#SF = float(SF)
#CN = input("Ratio C/N in Feed 1: ")
#CN = float(CN)
#D = input("Diluation rate of stage 1: ")
#D = float(D)

NF = SF*6*12/180 * CN

#Non-linear Equations (specific growth rate)
# m = (mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX))
# qC = (1.88*(1-rL)*bLC)
# bLC = (bCmax*1/(1+N/KiN)*(O/(KO2+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)*(1-C/Xf/KiC))
# bC = ((1-rL)*bLC)
# qL = (aL*m + bL)
# bL = (rL*bLC - KSL*(L/Xf+L)*(O/(KO2+O)))
# bE = (rE*bLC*(O/(KO2+O)) - KSE*(E/Xf+L)*(O/(KO2+O)))
# qE = (aE*m + bE)
# qS = (m/YXS + (O/(KO1*Xf+O))*(S/(KS+S))*mS + bC/YCS + (aL*m+bL)/YLS)
# qN = (m/YXN)

#Fb cal: eq1 = (V/1000) * ((7.14/YXN)*(mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX))*Xf + 1.59*((1-rL)*(bCmax*1/(1+N/KiN)*(O/(KO2+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)*(1-C/Xf/KiC)))*Xf) - FB

#Solve F, FS, FB
import numpy as np
from sympy import symbols, Eq, solve
F, FS, FB = symbols('F FS FB')
eq1 = Eq(FB*8/100, FS)
eq2 = Eq((FB+FS)/V, D)
eq3 = Eq(FS + FB, F)
sol = solve((eq1, eq2, eq3), (F, FS, FB))
F = sol[F]
FS = sol[FS]
FB = sol[FB]

#Solve Non-linear Equations
from scipy import stats
from scipy import optimize
def equations(vars):
    S, N, Xf, L, E, C  = vars
    eq4 = Xf*(-F+(mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX))*V)
    eq5 = FS*SF - F*S - ((mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX))/YXS + (O/(KO1*Xf+O))*(S/(KS+S))*mS + ((1-rL)*(bCmax*1/(1+N/KiN)*(O/(KO2+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)*(1-C/Xf/KiC)))/YCS + (aL*(mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX))+(rL*(bCmax*1/(1+N/KiN)*(O/(KO2+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)*(1-C/Xf/KiC)) - KSL*(L/Xf+L)*(O/(KO2+O))))/YLS)*Xf*V
    eq6 = FS*NF - F*N - ((mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX))/YXN)*Xf*V
    eq7 = -F*L + (aL*(mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)) + (rL*(bCmax*1/(1+N/KiN)*(O/(KO2+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)*(1-C/Xf/KiC)) - KSL*(L/Xf+L)*(O/(KO2+O))))*Xf*V
    eq8 = -F*E + (aE*(mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)) + (rE*(bCmax*1/(1+N/KiN)*(O/(KO2+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)*(1-C/Xf/KiC))*(O/(KO2+O)) - KSE*(E/Xf+L)*(O/(KO2+O))))*Xf*V
    eq9 = L - C
#eq9 = -F*C + (1.88*(1-rL)*(bCmax*1/(1+N/KiN)*(O/(KO2+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)*(1-C/Xf/KiC)))*Xf*V
    return (eq4, eq5, eq6, eq7, eq8, eq9)
def f(vars):
    return abs(np.array(sum(equations(vars))**2)-0)

S, N, Xf, L, E, C =  optimize.fsolve(equations, (optimize.fmin(f, (3.3, 0.05, 100, 70, 30, 60))))

#Result Validation
print(Xf*(-F+(mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX))*V))
print(FS*SF - F*S - ((mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX))/YXS + (O/(KO1*Xf+O))*(S/(KS+S))*mS + ((1-rL)*(bCmax*1/(1+N/KiN)*(O/(KO2+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)*(1-C/Xf/KiC)))/YCS + (aL*(mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX))+(rL*(bCmax*1/(1+N/KiN)*(O/(KO2+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)*(1-C/Xf/KiC)) - KSL*(L/Xf+L)*(O/(KO2+O))))/YLS)*Xf*V)
print(FS*NF - F*N - ((mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX))/YXN)*Xf*V)
print(-F*L + (aL*(mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)) + (rL*(bCmax*1/(1+N/KiN)*(O/(KO2+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)*(1-C/Xf/KiC)) - KSL*(L/Xf+L)*(O/(KO2+O))))*Xf*V)
print(-F*E + (aE*(mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)) + (rE*(bCmax*1/(1+N/KiN)*(O/(KO2+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)*(1-C/Xf/KiC))*(O/(KO2+O)) - KSE*(E/Xf+L)*(O/(KO2+O))))*Xf*V)
print(L-C)


m = (mmax*(N/(KN+N))*(O/(KO1+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX))
bLC = (bCmax*1/(1+N/KiN)*(O/(KO2+O))*(S/(KS+S))*1/(1+S/KiS)*1/(1+Xf/KiX)*(1-C/Xf/KiC))
qC = (1.88*(1-rL)*bLC)
bC = ((1-rL)*bLC)
bL = (rL*bLC - KSL*(L/Xf+L)*(O/(KO2+O)))
qL = (aL*m + bL)
bE = (rE*bLC*(O/(KO2+O)) - KSE*(E/Xf+L)*(O/(KO2+O)))
qE = (aE*m + bE)
qS = (m/YXS + (O/(KO1*Xf+O))*(S/(KS+S))*mS + bC/YCS + (aL*m+bL)/YLS)
qN = (m/YXN)

print("S = ", S, " N = ", N, " Xf = ", Xf," L = ", L, " E = ", E, " C = ", C,  )
     