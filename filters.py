"""
Design Filters for performing a rotationally invariant measurement of the
derivative.

This should be done following the procedure outlined in:

"""

def x_deriv(image, n_taps = 5):
    """
    Calculate the x-derivative of an image by convolving with the filter
    described at:

    Parameters:
    -----------
    image : np array of floating point values; shape: [n_rows, n_cols]

    Returns:
    ---------
    The x-derivative of the image
    """
    if n_taps = 5:
         k = np.array([0.030320, 0.249724, 0.439911, 0.249724, 0.030320])
         d = np.array([0.104550, 0.292315, 0.000000, -0.292315, -0.104550])
         kernel = np.outer(d, k)
    else:
         k = np.array([0.004711, 0.069321, 0.245410, 0.361117,
                       0.245410, 0.069321, 0.004711])
         d = np.array([0.018708, 0.125376, 0.193091, 0.000000, -0.193091,
                       -0.125376, -0.018708])
         kernel = np.outer(d, k)

def y_deriv(image, n_taps=5):
    """
    Calculate the x-derivative of an image by convolving with the filter
    described at:

    Parameters:
    -----------
    image : np array of floating point values; shape: [n_rows, n_cols]

    Returns:
    ---------
    The y-derivative of the image
    """
    if n_taps = 5:
         k = np.array([0.030320, 0.249724, 0.439911, 0.249724, 0.030320])
         d = np.array([0.104550, 0.292315, 0.000000, -0.292315, -0.104550])
         kernel = np.outer(k, d)
    else:
         k = np.array([0.004711, 0.069321, 0.245410, 0.361117,
                       0.245410, 0.069321, 0.004711])
         d = np.array([0.018708, 0.125376, 0.193091, 0.000000, -0.193091,
                       -0.125376, -0.018708])
         kernel = np.outer(k, d)

    return ndi.convolve(image, kernel, mode='constant')

def gabor_deriv(image, frequency, theta, sigma_x, sigma_y):
    """
    Convolves image with gabor filter that is specified by the parameters

    Parameters:
    -----------
    image : np array of floating point values; shape: [n_rows, n_cols]

    Returns:
    --------
    """
    kernel = gabor_kernel(frequency, theta, sigma_x, sigma_y)
    return ndi.convolve(image, kernel, mode='constant')
