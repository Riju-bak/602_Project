# Kinetic Parameters
KiC = 1.15
KiN = 0.073
KiX = 152.0
KiS = 100.0
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
O = 25.0
SF = 700.0
CN = 68.0
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
    S, N, Xf, L, E, C = vars
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
    eq4 = -F*Xf + m*Xf*V
    eq5 = FS*SF - F*S - qS*Xf*V
    eq6 = FS*NF - F*N - qN*Xf*V
    eq7 = -F*L + qL*Xf*V
    eq8 = -F*E + qE*Xf*V
    # eq9 = -F*C + qC*Xf*V
    eq9 = L - C
    return [eq4, eq5, eq6, eq7, eq8, eq9]

def f(vars):
    print("############ ITERATION ##############")
    print("xk: ", vars)
    print("fk: ", sum(np.array(equations(vars)) ** 2))
    print("####################################")
    return sum(np.array(equations(vars)) ** 2)

f_op_one = optimize.fmin(f, (3.3, 0.05, 100.0, 70.0, 30.0, 60.0), xtol=1e-4, maxiter=10000)
S, N, Xf, L, E, C =  optimize.fsolve(equations, f_op_one)
# S, N, Xf, L, E, C = f_op_one

print("S = ", S, " N = ", N, " Xf = ", Xf," L = ", L, " E = ", E, " C = ", C,  )

# EPA_titer = E/X*100
EFT = V/F # Effective fermentation time(h)
EPA_prod = E/EFT

X = Xf + L
EPA_percent_Biomass = E * 100 / X
Lipid_percent_Biomass = L*100/X
EPA_percent_Lipid = E*100/L

print("EPA_percent_Biomass: ", EPA_percent_Biomass)
print("Lipid_percent_Biomass: ", Lipid_percent_Biomass)
print("EPA_percent_Lipid: ", EPA_percent_Lipid)
print("EPA_productivity: ", EPA_prod)