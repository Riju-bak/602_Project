import numpy as np
from scipy import stats
from scipy import optimize
from pickle import TRUE
import pandas as pd

def continuosSolver():
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


    # Operating Conditions Input
    n = int(input("Number of stage: "))
    O = 25

    # n-stage
    data = []
    result = []
    Scon = 0
    for a in range (0, n+1, 1):
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
            FB = 0
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
                print("Volume of stage {}-reactor = ". format(a), V)
            F = V*D
        name = ['F','FB', 'SF', 'D', 'V', 'S', 'N', 'Xf', 'L', 'E', 'C']
        value = [F, FB, SF, D, V, S, N, Xf, L, E, C]
        data.append((dict(zip(name, value))))
    df = pd.DataFrame(data)

    batch_yield = float(input("STD Batch conversion yield: "))

         #Solve Non-linear Equations
    for a in range (1, n+1, 1):
        def equations(vars):
            FB, S, N, Xf, L, E, C  = vars

            #Specific growth rates
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

            #Feeds
            F = df.loc[(a), 'F']
            V = df.loc[(a), 'V']
            D = df.loc[(a), 'D']
            FS = D*V - FB/1000 #FB(mL/h)

            #Equations
            eq1 = FB - (V) * ((7.14/YXN)*m*Xf + 1.59*bC*Xf)
            eq2 = df.loc[(a-1), 'F']*df.loc[(a-1), 'Xf'] - F*Xf + m*Xf*V
            eq3 = df.loc[(a-1), 'F']*df.loc[(a-1), 'S'] + FS*df.loc[(a), 'SF'] - F*S - qS*Xf*V
            if a == 1:
                eq4 = df.loc[(a-1), 'F']*df.loc[(a-1), 'N'] + FS*NF - F*N - qN*Xf*V
            else:
                eq4 =  df.loc[(a-1), 'F']*df.loc[(a-1), 'N'] - F*N - qN*Xf*V
            eq5 = df.loc[(a-1), 'F']*df.loc[(a-1), 'L'] - F*L + qL*Xf*V
            eq6 = df.loc[(a-1), 'F']*df.loc[(a-1), 'E'] - F*E + qE*Xf*V
            eq7 = df.loc[(a-1), 'F']*df.loc[(a-1), 'C'] - F*C + qC*Xf*V
            return (eq1, eq2, eq3, eq4, eq5, eq6, eq7)

        def f(vars):
            return abs(np.array(sum(equations(vars))**2)-0)

        #Initial Guesses
        IGuess1 = np.array([1, 0.05, 0.01, 150, 100, 50, 100])
        IGuess2 = np.array([1, 0.05, 0.001, 70, 100, 50, 100])
        if a in range (1, 3, 1):
            FB, S, N, Xf, L, E, C = optimize.fsolve(equations, (optimize.fmin(f, IGuess1, xtol=0.01, maxiter=10000)))
        else:
            FB, S, N, Xf, L, E, C = optimize.fsolve(equations, (optimize.fmin(f, IGuess2, xtol=0.01, maxiter=10000)))

        roots = FB, S, N, Xf, L, E, C
        df.loc[a, 'FB'], df.loc[a, 'S'], df.loc[a, 'N'], df.loc[a, 'Xf'], df.loc[a, 'L'], df.loc[a, 'E'], df.loc[a, 'C'] = FB, S, N, Xf, L, E, C

        #Validate results
        if a in range (1, 3, 1):
            print(optimize.fmin(f, IGuess1, xtol=0.01, maxiter=10000))
        else:
            print(optimize.fmin(f, IGuess2, xtol=0.01, maxiter=10000))
        print(equations(roots))
        print(" ")
        print(" ")

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

        F = df.loc[(a), 'F']
        V = df.loc[(a), 'V']
        D = df.loc[(a), 'D']
        FS = D*V - FB/1000

        #OUR
        qO2 = -(1.07*qS - 1.37*m - 2.86*(aL*m+bL) - 1.45*bC)
        OUR = 31.25*qO2*Xf
        print("OUR of stage ", a, "(mmol/h): ", OUR)

        #Evaluation results
        stage = a
        t = 1/D
        EPA_titer = E*100/(Xf+L)
        EPA_rate = E/t
        Scon = Scon + FS*df.loc[(a), 'SF'] - F*S
        EPA_in_Lipid = E*100/L
        Lipid_content = L*100/(Xf+L)
        A = ['Stage', 't','EPA Titer', 'EPA Rate', 'Scon','EPA in Lipid', 'Lipid Content' ]
        B = [stage, t, EPA_titer, EPA_rate, Scon, EPA_in_Lipid, Lipid_content]
        result.append((dict(zip(A, B))))
    results = pd.DataFrame(result)

    # Overall result
    t = df['V'].sum()/F
    EPA_titer = E*100/(Xf+L)

    EPA_rel_yield = (E*F/Scon)*(100/batch_yield)

    EPA_in_Lipid = E*100/L
    Lipid_content = L*100/(Xf+L)

    sum_EPA_rate = 0
    for a in range (1, n+1, 1):
        EPA_rate = df.loc[a, 'V']*results.loc[(a-1), 'EPA Rate']
        sum_EPA_rate = sum_EPA_rate + EPA_rate
    EPA_rate = sum_EPA_rate/df['V'].sum()

    overall = {'Stage': 'Overall', 't':t, 'EPA Titer':EPA_titer, 'EPA Rate':EPA_rate,'EPA in Lipid':EPA_in_Lipid, 'Lipid Content':Lipid_content}
    results = results.append(overall, ignore_index=True)


    df = np.round(df, decimals = 4)
    results = np.round(results, decimals = 2)
    print(" ")
    print(" ")
    print(" ")
    print("Species Concentration")
    print(df)
    #print(df[['D', 'V', 'S', 'N', 'Xf', 'L', 'E', 'C']])
    print(" ")
    print(" ")
    print(" ")
    print("Evaluation results")
    print(results)
    print(EPA_rel_yield)

