from utils import find_ind_val
import numpy as np
from scipy.integrate import simpson


def EPAYieldSolver(E, S, qSXf_integrated, t, **kwargs):
    """This method returns the amount of EPA produced per amount of S consumed (unit/g)"""
    S_20_ind = find_ind_val(S, 20)  # The index at which S becomes 20, used to find the time at which S becomes 20
    qSXf = np.diff(qSXf_integrated) / np.diff(t)
    qsXf_extra = qSXf[S_20_ind:]
    t_ = np.linspace(t[S_20_ind], t[-1], len(qsXf_extra))
    S_consumed_extra = simpson(qSXf[S_20_ind:], t_)  # Using simpson integration to find amount of S that reacts
    # after we start feeding glucose
    S_consumed = (S[0] - 20) + S_consumed_extra
    EPA_produced = (E[-1] - E[0])
    EPA_yield = EPA_produced / S_consumed
    return EPA_yield
