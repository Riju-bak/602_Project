import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from utils import plot_2D
from .epa_yield import EPAYieldSolver

np.seterr(divide='ignore', invalid='ignore')


def batchIVPSolver():
    # Kinetic parameters
    KiC = 1.15
    KiN = 0.073
    KiX = 152
    KiS = 100
    KN = 0.033
    KO1 = 0.65
    KO2 = 5.5
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
    bCmax = 0.09
    mu_max = 0.26

    # Init vals
    S0 = 47.6
    N0 = 7.58
    O0 = 25
    X0 = 1
    SF = 700
    V0 = 1.05
    L0 = 0
    E0 = 0
    LE0 = 0
    C0 = 0
    OUR0 = 0
    O = O0

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

        return [dXfdt, dCdt, dLdt, dLEdt, dEdt, dSdt, dNdt, dVdt, FS]

    # Initial Conditions
    # FS(t=0) = 0
    y0 = [X0, C0, L0, LE0, E0, S0, N0, V0, 0]

    # Test ODEs
    # print(odes(y=y0, t=0))

    # declare a time vector(time window)
    T = 144  # Total batch operation time(h)
    T_split = 10000
    t = np.linspace(0, T, T_split)
    t_span = (0, T)

    y = solve_ivp(odes, t_span, y0, t_eval=t, method='LSODA')

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
    EFT = 120  # Effective fermentation time(h)
    EPA_prod = (E - E[0]) / (t - t[0])

    L = LE + E
    X = Xf + L  # Biomass (Unit/L)

    EPA_percent_Biomass = E * 100 / X
    Lipid_percent_Biomass = L * 100 / X
    EPA_percent_Lipid = E * 100 / L

    # EPA produced per amount of S consumed (unit/g)
    FS_integrated = y.y[8, :]
    EPA_yield, S_consumed_extra, S_consumed, EPA_produced = EPAYieldSolver(E, S, FS_integrated, V, SF, t)
    print("EPA_yield: ", EPA_yield)
    ###################################################

    # Plotting X
    plot_2D(X, t, figure_name="Biomass", title="Biomass(Unit/L)", ylabel="X", xlabel="time(h)")

    # Plotting lipid(% Biomass)
    plot_2D(Lipid_percent_Biomass, t, figure_name="lipid(% Biomass)", title="lipid(% Biomass)", xlabel="time(h)")

    # Plotting EPA(% Liquid)
    plot_2D(EPA_percent_Lipid, t, figure_name="EPA(% Lipid)", title="EPA(% Lipid)", xlabel="time(h)")

    # Plotting EPA(% Biomass)
    plot_2D(EPA_percent_Biomass, t, figure_name="EPA(% Biomass)", title="EPA(% Biomass)", xlabel="time(h)")

    # Plotting EPA_prod(unit/L/h)
    plot_2D(EPA_prod, t, figure_name="EPA_prod", title="EPA_prod(Unit/L/h)", xlabel="time(h)")

    # Plotting S
    plot_2D(S, t, figure_name="S vs t", title="S vs t", ylabel="S", xlabel="time(h)")

    # Plotting N
    plot_2D(N, t, figure_name="N vs t", title="N vs t", ylabel="N", xlabel="time(h)")

    # Plotting E
    plot_2D(E, t, figure_name="EPA_conc", title="EPA_conc(Unit/L)", ylabel="E", xlabel="time(h)")

    plt.show()
