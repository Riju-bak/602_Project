from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from dotenv import load_dotenv
import os


def odes(y, t):
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

    # assign each ODE to a vector elem
    Xf = y[0]
    C = y[1]
    L = y[2]
    LE = y[3]
    E = y[4]
    S = y[5]
    N = y[6]
    V = y[7]

    ########## ODE Definitions #####################
    FB = V / 1000 * (7.14 / YXN * (mu_max * (N / (KN + N)) * (O / (KO1 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf))) * Xf + 1.59 * ((1 - rL) * (bCmax * (KiN / (KiN + N)) * (O / (KO2 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf)) * ((KiC - C / Xf) / KiC))) * Xf)

    bL = rL * (bCmax * (KiN / (KiN + N)) * (O / (KO2 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf)) * ((KiC - C / Xf) / KiC)) - KSL * (L / (Xf + L)) * (O / (KO2 + O))
    qS = (mu_max * (N / (KN + N)) * (O / (KO1 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf))) / YXS + (O / (KO1 * Xf + O)) * (S / (KS + S)) * mS + ((1 - rL) * (bCmax * (KiN / (KiN + N)) * (O / (KO2 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf)) * ((KiC - C / Xf) / KiC))) / YCS + (aL * (mu_max * (N / (KN + N)) * (O / (KO1 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf))) + bL) / YLS

    FS = 0 if S>20 else (qS * Xf + FB * S / V) * (V / (SF - S))  # Glucose Feed Rate (L/h)
    D = (FB + FS) / V

    dXfdt = (mu_max * (N / (KN + N)) * (O / (KO1 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf))) * Xf - D * Xf

    dCdt = (1.88 * (1 - rL) * (bCmax * (KiN / (KiN + N)) * (O / (KO2 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf)) * ((KiC - C / Xf) / KiC))) * Xf - D * C

    qL = aL * (mu_max * (N / (KN + N)) * (O / (KO1 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf))) + bL
    dLdt = qL * Xf - D * L
    bE = rE * (bCmax * (KiN / (KiN + N)) * (O / (KO2 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf)) * ((KiC - C / Xf) / KiC)) * (O / (KO2 + O)) - KSE * (E / (Xf + L)) * (O / (KO2 + O))
    dLEdt = ((aL - aE) * (mu_max * (N / (KN + N)) * (O / (KO1 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf))) + (bL - bE)) * Xf - D * LE

    qE = aE * (mu_max * (N / (KN + N)) * (O / (KO1 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf))) + bE
    dEdt = qE * Xf - D * E

    dSdt = (-qS * Xf - (FB * S/V)) if S>20 else 0

    qN = (mu_max * (N / (KN + N)) * (O / (KO1 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf))) / YXN
    dNdt = qN * Xf + D * N

    dVdt = FS + FB

    qO2 = -(1.07 * qS - 1.37 * (mu_max * (N / (KN + N)) * (O / (KO1 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf))) - 2.86 * (aL * (mu_max * (N / (KN + N)) * (O / (KO1 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf))) + bL) - 1.45 * ((1 - rL) * (bCmax * (KiN / (KiN + N)) * (O / (KO2 + O)) * (S / (KS + S)) * (KiS / (KiS + S)) * (KiX / (KiX + Xf)) * ((KiC - C / Xf) / KiC))))
    OUR = 31.25 * qO2 * Xf

    return [dXfdt, dCdt, dLdt, dLEdt, dEdt, dSdt, dNdt, dVdt]
###############################################################################################################################################

load_dotenv()

S0 = float(os.environ.get("S0"))
N0 = float(os.environ.get("N0"))
X0 = float(os.environ.get("X0"))
V0 = float(os.environ.get("V0"))
L0 = float(os.environ.get("L0"))
E0 = float(os.environ.get("E0"))
LE0 = float(os.environ.get("LE0"))
C0 = float(os.environ.get("C0"))

# Initial Conditions
y0 = [X0, C0, L0, LE0, E0, S0, N0, V0]

# Test ODEs
# print(odes(y=y0, t=0))

# declare a time vector(time window)
t = np.linspace(0, 144, 1000)
y = odeint(odes, y0, t)

Xf = y[:, 0]
C = y[:, 1]
L = y[:, 2]
LE = y[:, 3]
E = y[:, 4]
S = y[:, 5]
N = y[:, 6]
V = y[:, 7]

# EPA_titer = E/X*100
# EPA_prod = E/t where t=effective fermentation time
L = LE + E
X = Xf + L

EPA_percent_Biomass = E * 100 / X
Lipid_percent_Biomass = L*100/X

#Plotting X
plt.figure()
plt.title("Biomass(Unit/L)")
plt.ylabel("X")
plt.xlabel("t")
plt.plot(t, X)
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

#Plotting EPA(% Biomass)
plt.figure()
plt.title("EPA(% Biomass)")
plt.xlabel("time(h)")
plt.plot(t, EPA_percent_Biomass)
plt.draw()

#Plotting lipid(% Biomass)
plt.figure()
plt.title("lipid(% Biomass)")
plt.xlabel("time(h)")
plt.plot(t, Lipid_percent_Biomass)
plt.draw()

plt.show()


