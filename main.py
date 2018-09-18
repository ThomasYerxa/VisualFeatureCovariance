"""Main loop for measuring orientation"""

import numpy as np
from utils import (load_data, lms_to_xyz, normalize, gen_pyramid, gen_cov,
                   extract_features)

ORIENTEDNESS_THRESHOLD = 0.8 #Single value threshold
ENERGY_THRESHOLD = 0.68 #Percentile threshold
n_scales = 6 #number of spatial scales to examine

#Load data, convert to XYZ, isolate luminance channel
RGB = load_data()
Y = []
for i in range(len(RGB)):
    X = rgb_to_xyz(RGB[i])
    Y.append(X[:,:,1])
Y = np.array(Y)

n_images, n_rows, n_cols = Y.shape

def calc_orientation_features(image, n_taps=5):
    """
    Calculates the orientation tensor features for each position in an image

    Parameters:
    -----------
    image : np array of floating point values; shape: [n_rows, n_cols]
    n_taps : specifies whether to use 5-tap or 7-tap derivative filters

    Returns:
    --------
    E : list of energy values for each pixel position in image
    O : list of orientedness values for each pixel position in image
    T : list of theta values for each pixel position in image
    """
    E = [] #energy values
    O = [] #orientedness values
    T = [] #theta values

    n_rows, n_cols = iamge.shape
    for x in range(n_cols):
        for y in range(n_rows):
            x_p = x_deriv(image, n_taps=n_taps)
            y_p = y_deriv(image, n_taps=n_taps)
            T_loc = gen_cov(x_p, y_p, x, y)
            energy, orientedness, theta = extract_features(T_loc)
            E.append(energy)
            O.append(orientedness)
            T.append(theta)

    return E, O, T


def calc_gabor_features(frequencies, orientations, stds):
    """
    *IN PROGRESS*  (Not sure about approriate feature set given filtered image)

    Calculates the orientation distribution using the gabor filter approach
    """
    means = []
    vars = []
    for i in range(n_images):
        for f in frequencies:
            for theta in orientations:
                for std in stds:
                    filtered = gabor_deriv(Y[i], f, t, std, std)
                    means.append(np.mean(filtered))
                    vars.append(np.var(filtered))

    return maens, vars

def main_one(scale):
    """
    *IN PROGRESS*

    Calculates the orientation distribution using the rotationally invariant
    derivative approach
    """
    E = [] #energy values
    O = [] #orientedness values
    T = [] #theta values
    for i in range(Y.shape[0]):
        image = gen_pyramid(Y[i], scale)[-1]
        E_i, O_i, T_i = calc_orientation_features(image)
        E += E_i
        O += O_i
        T += T_i

    E = np.array(E)
    O = np.array(O)
    T = np.array(T)

    E_thresh = np.percentile(E, 0.68)
    theta_significant = []

    for i in range(E.shape[0]):
        if O[i] >= 0.8 & E[i] > E_thresh:
            theta_significant.append(T[i])


def main_two():
    """
    *NOT IMPLEMENTED*
    """
    return None
