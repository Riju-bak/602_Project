from .batch_solve_ivp import batchIVPSolver
import pandas as pd


def generate_bsivp_report():
    f = open("./bsivp_report.txt", "w")
    solverMethods = {
        0: "RK45",
        1: "RK23",
        2: "DOP853",
        3: "Radau",
        4: "BDF",
        5: "LSODA",
    }
    dct = {"Method": {}, "EPA Yield": {}}
    for solver_id in solverMethods.keys():
        EPA_yield = batchIVPSolver(generate_report=True, solver_id=solver_id)
        dct["Method"][solver_id] = solverMethods[solver_id]
        dct["EPA Yield"][solver_id] = EPA_yield
    data = pd.DataFrame(dct)
    print(data)
    f.write(str(data))
    f.close()