import numpy as np
from .f90wrap.qcd import dummy as f90src

def alphaQCD(Q2):
    if(np.isscalar(Q2)):
        Q2 = np.array([Q2])
    return f90src.alpha_wrap(Q2)
