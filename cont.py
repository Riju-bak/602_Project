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

SF = float(input("Glucose Concentration in Feed 1: "))
CN = float(input("Ratio C/N in Feed 1: "))
D = float(input("Diluation rate of stage 1: "))

NF = (SF*6*12/180) / CN

#Solve Non-linear Equations - Stage 1
import numpy as np
from scipy import stats
from scipy import optimize
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
    eq1 = FB - ((V/1000) * ((7.14/YXN)*m*Xf + 1.59*bC*Xf))
    FS = D*V - FB
    F = D*V
    eq2 = -F*Xf + m*Xf*V
    eq3 = FS*SF - F*S - qS*Xf*V
    eq4 = FS*NF - F*N - qN*Xf*V
    eq5 = -F*L + qL*Xf*V
    eq6 = -F*E + qE*Xf*V
    eq7 = -F*C + qC*Xf*V
    return (eq1, eq2, eq3, eq4, eq5, eq6, eq7)
def f(vars):
    return abs(np.array(sum(equations(vars))**2)-0)
FB, S, N, Xf, L, E, C =  optimize.fsolve(equations, (optimize.fmin(f, (0.02, 0.095, 0.01, 150, 100, 56, 100), xtol=0.01, maxiter=100)))
print("FB = ", FB*1000, " FS = ", (D*V-FB)*1000," S = ", S, " N = ", N, " Xf = ", Xf," L = ", L, " E = ", E, " C = ", C, " t = ", 1/D )

#n-stage


#Validate results 
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
FS = D*V - FB
F = D*V


print(-F*Xf + m*Xf*V)
print(FS*SF - F*S - qS*Xf*V)
print(FS*NF - F*N - qN*Xf*V)
print(-F*L + qL*Xf*V)
print(-F*E + qE*Xf*V)
print(-F*C + qC*Xf*V)

EFT = V/F # Effective fermentation time(h)
EPA_rate = E/EFT

#Evaluation Results 
X = Xf + L
EPA_titer = E * 100 / X
Lipid_percent_Biomass = L*100/X
EPA_percent_Lipid = E*100/L

print("EPA titer: ", EPA_titer)
print("Lipid_percent_Biomass: ", Lipid_percent_Biomass)
print("EPA_percent_Lipid: ", EPA_percent_Lipid)
print("EPA_rate: ", EPA_rate)
