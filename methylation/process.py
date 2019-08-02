import numpy as np
import pandas as pd
import re


def m_from_beta(beta):
    return np.log2(beta / (1 - beta))


def beta_from_m(m):
    return 2. ** m / (1 + 2. ** m)
