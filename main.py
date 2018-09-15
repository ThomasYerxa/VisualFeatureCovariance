"""Main loop for measuring orientation"""

import numpy as np
from utils import (load_data, lms_to_xyz, normalize, gen_pyramid, gen_cov,
                   extract_features)

ORIENTEDNESS_THRESHOLD = 0.8 #Single value threshold
ENERGY_THRESHOLD = 0.68 #Percentile threshold
n_scales = 6 #number of spatial scales to examine

#Load data, convert to XYZ, isolate luminance channel
XYZ = load_data()
Y = []
for i in range(len(X)):
    X[i] = lms_to_xyz(X[i])
    Y.append(X[i][:,:,1])
Y = np.array(Y)


def gen_hist(x):
    """

    """
    for i in range(len(x)):


def main():
    """
    """
    scale = 0
