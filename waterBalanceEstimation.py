import pandas as pd
import numpy as np


def waterBalance(E, I, P, dV, deltaRin, deltaI, deltaP, deltaQ, deltaE):
    """
    Parameters:
        - E, I, P, dV in mm
        - the corresponding isotopic composition in per mill
    """
    Rin = (E*(deltaE-deltaQ) - I*(deltaI-deltaQ) - P*(deltaP-deltaQ))/(deltaRin-deltaQ)
    x = E/(I+P+Rin)
    Rout = (1-x)*(Rin+I+P) - Q - dV

    return Rin, Rout




