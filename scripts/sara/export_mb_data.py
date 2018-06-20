import pandas as pd
from rnaseq import loader
import os


if __name__ == "__main__":
    obj_star = loader.load_by_patient(['ICb1299', '3021'], source='star', type='cell_culture')
    obj_salmon = loader.load_by_patient(['ICb1299', '3021'], source='salmon', type='cell_culture')