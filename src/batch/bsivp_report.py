from .batch_solve_ivp import batchIVPSolver
import pandas as pd


def generate_yield_report():
    f = open("./yield_report.txt", "w+")
    solverMethods = {
        0: "RK45",
        1: "RK23",
        2: "DOP853",
        3: "Radau",
        4: "BDF",
        5: "LSODA",
    }
    dct = {"Method": {}, "EPA_Yield": {}, "S_consumed_extra": {}, "EPA_produced": {} }
    for solver_id in solverMethods.keys():
        EPA_yield, S_consumed_extra, EPA_produced = batchIVPSolver(generate_report=True, solver_id=solver_id)
        dct["Method"][solver_id] = solverMethods[solver_id]
        dct["EPA_Yield"][solver_id] = EPA_yield
        dct["S_consumed_extra"][solver_id] = S_consumed_extra
        dct["EPA_produced"][solver_id] = EPA_produced
    data = pd.DataFrame(dct)
    print(data)
    f.write(str(data))
    f.close()

def generate_batch_results_report():
    f = open("./batch_results.txt", "w+")
    # EPA_titer, EPA_rate, EPA_Yield = batchIVPSolver(generate_batch_result_report=True)
    solverMethods = {
        0: "RK45",
        1: "RK23",
        2: "DOP853",
        3: "Radau",
        4: "BDF",
        5: "LSODA",
    }
    dct = {"Method": {}, "EFT(h)": {}, "Glucose-Feed(g/L)": {}, "EPA_Titer": {}, "EPA_Rate": {}, "EPA_Yield": {},
           "EPA_conc": {}, "Biomass_conc": {}, "Lipid_Titer": {}, "EPA_percent_Lipid": {}}
    for solver_id in solverMethods.keys():
        results = batchIVPSolver(generate_batch_result_report=True, solver_id=solver_id)
        dct["Method"][solver_id] = solverMethods[solver_id]
        dct["EFT(h)"][solver_id] = 120
        dct["Glucose-Feed(g/L)"][solver_id] = 700
        for k in results.keys():
            dct[k][solver_id] = results[k]
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 150)
    data = pd.DataFrame(dct)
    print(data)
    f.write(str(data))
    f.close()

