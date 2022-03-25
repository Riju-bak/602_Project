import os
import matplotlib.pyplot as plt
import numpy as np
from dotenv import load_dotenv
from scipy.integrate import solve_ivp
from utils import find_ind_val

np.seterr(divide='ignore', invalid='ignore')

load_dotenv()

# Kinetic Params
KiC = float(os.environ.get("KiC"))
KiN = float(os.environ.get("KiN"))
KiX = float(os.environ.get("KiX"))
KiS = float(os.environ.get("KiS"))
KN = float(os.environ.get("KN"))
KO1 = float(os.environ.get("KO1"))
KO2 = float(os.environ.get("KO2"))
KS = float(os.environ.get("KS"))
KSE = float(os.environ.get("KSE"))
KSL = float(os.environ.get("KSL"))
mS = float(os.environ.get("mS"))
rE = float(os.environ.get("rE"))
rL = float(os.environ.get("rL"))
YCS = float(os.environ.get("YCS"))
YLS = float(os.environ.get("YLS"))
YXN = float(os.environ.get("YXN"))
YXS = float(os.environ.get("YXS"))
aE = float(os.environ.get("aE"))
aL = float(os.environ.get("aL"))
bCmax = float(os.environ.get("bCmax"))
mu_max = float(os.environ.get("mu_max"))

# Init vals
SF = float(os.environ.get("SF"))
O0 = float(os.environ.get("O0"))
O = O0  # O is assumed constant

S0 = float(os.environ.get("S0"))
N0 = float(os.environ.get("N0"))
X0 = float(os.environ.get("X0"))
V0 = float(os.environ.get("V0"))
L0 = float(os.environ.get("L0"))
E0 = float(os.environ.get("E0"))
LE0 = float(os.environ.get("LE0"))
C0 = float(os.environ.get("C0"))


# OUR0 = float(os.environ.get("OUR0"))

def odes(t, y):
    # assign each ODE to a vector elem
    Xf = y[0]
    C = y[1]
    L = y[2]
    LE = y[3]
    E = y[4]
    S = y[5]
    N = y[6]
    V = y[7]

    # ########## ODE Definitions #####################
    mu = (mu_max * (N / (KN + N)) * (O / (KO1 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf)))
    bLC = (bCmax * (KiN / (KiN + N)) * (O / (KO2 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf)) * (
                (KiC - C / Xf) / KiC))
    bC = ((1 - rL) * bLC)
    FB = V / 1000 * (7.14 / YXN * mu * Xf + 1.59 * bC * Xf)

    bL = rL * bLC - KSL * (L / (Xf + L)) * (O / (KO2 + O))
    qS = mu / YXS + (O / (KO1 * Xf + O)) * (S / (KS + S)) * mS + bC / YCS + (aL * mu + bL) / YLS

    FS = 0 if S>20 else (qS * Xf + FB * S / V) * (V / (SF - S))  # Glucose Feed Rate (L/h)
    D = (FB + FS) / V

    dXfdt = mu * Xf - D * Xf

    qC = (1.88 * (1 - rL) * bLC)
    dCdt = qC * Xf - D * C

    qL = aL * mu + bL
    dLdt = qL * Xf - D * L
    bE = rE * bLC * (O / (KO2 + O)) - KSE * (E / (Xf + L)) * (O / (KO2 + O))
    dLEdt = ((aL - aE) * mu + (bL - bE)) * Xf - D * LE

    qE = aE * mu + bE
    dEdt = qE * Xf - D * E

    dSdt = (-qS * Xf - (FB * S/V)) if S>20 else 0

    qN = mu / YXN
    dNdt = -qN * Xf - D * N

    dVdt = FS + FB

    ###############################################################################################################################################

    return [dXfdt, dCdt, dLdt, dLEdt, dEdt, dSdt, dNdt, dVdt]


# Initial Conditions
y0 = [X0, C0, L0, LE0, E0, S0, N0, V0]

# Test ODEs
# print(odes(y=y0, t=0))

# declare a time vector(time window)
T = 144  # Total batch operation time(h)
T_split = 10000
t = np.linspace(0, T, T_split)
t_span = (0, T)
solver_methods = ['RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA']
print("Choose Solver method: \n0)RK45\n1)RK23\n2)DOP853\n3)Radau\n4)BDF\n5)LSODA")
solve_ivp_method = solver_methods[int(input("Choice: "))]
y = solve_ivp(odes, t_span, y0, t_eval=t, method=solve_ivp_method)

Xf = y.y[0, :]
C = y.y[1, :]
L = y.y[2, :]
LE = y.y[3, :]
E = y.y[4, :]
S = y.y[5, :]
N = y.y[6, :]
V = y.y[7, :]
# OUR = y[:, 8]

# EPA_titer = E/X*100
EFT = 120 # Effective fermentation time(h)
EPA_prod = (E - E[0])/(t - t[0])

L = LE + E
X = Xf + L #Biomass (Unit/L)

EPA_percent_Biomass = E * 100 / X
Lipid_percent_Biomass = L*100/X
EPA_percent_Lipid = E*100/L

#EPA produced per amount of S consumed (unit/g)
qSXf_integrated = y[:, 8]
S_20_ind = find_ind_val(S, 20)
qSXf = np.diff(qSXf_integrated)/np.diff(t)
qsXf_extra = qSXf[S_20_ind:]
t_ = np.linspace(t[S_20_ind], T, len(qsXf_extra))
S_consumed_extra = simpson(qSXf[S_20_ind:], t_)
S_consumed = (S[0]-20) + S_consumed_extra
EPA_produced = (E[-1]-E[0])
EPA_yield = EPA_produced/S_consumed
print("EPA_yield: ", EPA_yield)
###################################################

#Plotting X
plt.figure()
plt.title("Biomass(Unit/L)")
plt.ylabel("X")
plt.xlabel("t")
plt.plot(t, X)
plt.draw()

#Plotting lipid(% Biomass)
plt.figure()
plt.title("lipid(% Biomass)")
plt.xlabel("time(h)")
plt.plot(t, Lipid_percent_Biomass)
plt.draw()

#Plotting OUR
# plt.figure()
# plt.title("OUR(mmol/L/h)")
# plt.xlabel("time(h)")
# plt.plot(t, OUR)
# plt.draw()

#Plotting EPA(% Liquid)
plt.figure()
plt.title("EPA(% Liquid)")
plt.xlabel("time(h)")
plt.plot(t, EPA_percent_Lipid)
plt.draw()

#Plotting EPA(% Biomass)
plt.figure()
plt.title("EPA(% Biomass)")
plt.xlabel("time(h)")
plt.plot(t, EPA_percent_Biomass)
plt.draw()

#Plotting EPA_prod(unit/L/h)
plt.figure()
plt.title("EPA_productivity (Unit/L/h)")
plt.xlabel("time(h)")
plt.plot(t, EPA_prod)
plt.draw()

#Plotting S
plt.figure()
plt.title("S vs t")
plt.ylabel("S")
plt.xlabel("t")
plt.plot(t, S)
plt.draw()

#Plotting N
plt.figure()
plt.title("N vs t")
plt.ylabel("N")
plt.xlabel("t")
plt.plot(t, N)
plt.draw()

# Plotting E
plt.figure()
plt.title("E vs t")
plt.ylabel("E")
plt.xlabel("t")
plt.plot(t, E)
plt.draw()


plt.show()
