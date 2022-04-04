from utils import find_ind_val
import numpy as np
from scipy.integrate import simpson, romb

def EPAYieldSolver(E, S, FS_integrated, V, SF, t, **kwargs):
    """This method returns the amount of EPA produced per amount of S consumed (unit/g)"""
    S_20_ind = find_ind_val(S, 20)  # The index at which S becomes 20, used to find the time at which S becomes 20
    FS = np.diff(FS_integrated) / np.diff(t)
    FS = FS[S_20_ind:]
    t_ = np.linspace(t[S_20_ind], t[-1], len(FS))
    S_consumed_extra = SF*simpson(FS, t_)  # Using simpson integration to find amount of S that reacts
    # # after we start feeding glucose
    S_consumed = (S[0]*V[0] - S[-1]*V[-1]) + S_consumed_extra
    EPA_produced = E[-1]*V[-1]
    EPA_yield = EPA_produced / S_consumed
    return EPA_yield, S_consumed_extra, S_consumed, EPA_produced