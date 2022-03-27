from scipy.integrate import odeint, simpson
import numpy as np
import matplotlib.pyplot as plt
import os
from dotenv import load_dotenv
from utils import find_ind_val, plot_2D
from .epa_yield import EPAYieldSolver

np.seterr(divide='ignore', invalid='ignore')


def batchSolver():
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

    def odes(y, t):
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

        FS = 0 if S > 20 else (qS * Xf + FB * S / V) * (V / (SF - S))  # Glucose Feed Rate (L/h)
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

        dSdt = (-qS * Xf - (FB * S / V)) if S > 20 else 0

        qN = mu / YXN
        dNdt = -qN * Xf - D * N

        dVdt = FS + FB

        ###############################################################################################################################################

        return [dXfdt, dCdt, dLdt, dLEdt, dEdt, dSdt, dNdt, dVdt, qS * Xf]

    # Initial Conditions
    # NOTE: qSXf(t=0) = 0
    y0 = [X0, C0, L0, LE0, E0, S0, N0, V0, 0]

    # Test ODEs
    # print(odes(y=y0, t=0))

    # declare a time vector(time window)
    T = 144  # Total batch operation time(h)
    T_split = 100000
    t = np.linspace(0, T, T_split)
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
    EFT = 120  # Effective fermentation time(h)
    EPA_prod = (E - E[0]) / (t - t[0])

    L = LE + E
    X = Xf + L  # Biomass (Unit/L)

    EPA_percent_Biomass = E * 100 / X
    Lipid_percent_Biomass = L * 100 / X
    EPA_percent_Lipid = E * 100 / L

    # EPA produced per amount of S consumed (unit/g)
    qSXf_integrated = y[:, 8]
    print("EPA_yield: ", EPAYieldSolver(E, S, qSXf_integrated, t))
    ###################################################

    # Plotting X
    plot_2D(X, t, figure_name="Biomass(Unit/L)", title="Biomass(Unit/L)", ylabel="X", xlabel="t")

    # Plotting lipid(% Biomass)
    plot_2D(Lipid_percent_Biomass, t, figure_name="lipid(% Biomass)", title="lipid(% Biomass)", xlabel="time(h)")

    # Plotting EPA(% Liquid)
    plot_2D(EPA_percent_Lipid, t, figure_name="EPA(% Liquid)", title="EPA(% Liquid)", xlabel="time(h)")

    # Plotting EPA(% Biomass)
    plot_2D(EPA_percent_Biomass, t, figure_name="EPA(% Biomass)", title="EPA(% Biomass)", xlabel="time(h)")

    # Plotting EPA_prod(unit/L/h)
    plot_2D(EPA_prod, t, figure_name="EPA_prod(Unit/L/h)", title="EPA_prod(Unit/L/h)", xlabel="time(h)")

    # Plotting S
    plot_2D(S, t, figure_name="S vs t", title="S vs t", ylabel="S", xlabel="t")

    # Plotting N
    plot_2D(N, t, figure_name="N vs t", title="N vs t", ylabel="N", xlabel="t")

    # Plotting E
    plot_2D(E, t, figure_name="E vs t", title="E vs t", ylabel="E", xlabel="t")

    plt.show()
