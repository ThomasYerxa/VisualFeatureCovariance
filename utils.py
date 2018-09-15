"""
The Goal of this program is to measure the distribution of orientation in
natural scenes by replicating the procedure outlined in:

These functions assume that 'images,' are in the LMS color space. For images
downloaded from:
in RGB, there is a matlab routine available to convert from RGB to LMS using
camera parameters specific to the data set.
"""

from sklearn.decomposition import PCA
import numpy as np

def load_data():
    """
    File I/O. Loads data (which should be images from LINK, that have been
    converted from RGB to LMS using the provided matlab routine)

    Returns:
    --------
    data : numpy array of floating point. shape: (n_images, im, im, 3)

    """
    return None

def lms_to_xyz(image):
    """
    Takes an LMS image to the XYZ space via a matrix multiplication.
    Parameters:
    -----------
    image : numpy array of floating point values. shape: (im, im, 3)

    Returns:
    --------
    im_new : numpy array of floating point values. shape: (im, im, 3)
    """

    return None

def normalize(image):
    """
    Normalizes an image by the mean luminance.

    Parameters:
    ------------
    image : numpy array of floating point values (corresponding to luminance
            data). shape: (im, im)

    Returns:
    ---------
    im_new : numpy array of floating point values (corresponding to normalized
             lumincance data). shape: (im, im)
    """

    return image / np.mean(image, axis=0)


def gen_pyramid(image, n_layers):
    """
    Generates a Gaussian Pyramid of images from normalized luminance data.
    The lowest layer is the image, successive layers capture spatial features
    of larger scales at lower resolutionself.

    Parameters:
    -----------
    image : numpy array of floating point values (corresponding to normalized
            lumincance data). shape: (im, im)

    Returns:
    --------
    pyr : a list of diff
    """
    # Discrete approximation to Gaussian Filter with sigma = 1
    G = [[1, 4, 7, 4, 1],
         [4, 16, 26, 16, 4],
         [7, 26, 41, 26, 41],
         [4, 16, 26, 16, 4],
         [1, 4, 7, 4, 1]
        ]
    G = np.array(G) * 1.0 / (273.0 )
    pyr = [image]
    for i in range(n_layers):
        # convolve pyr[-1] with gaussian filter, then downsample the blurred
        # image to 1/2 the original size
        smoothed = scipy.signal.convolve2d(pyr[-1], G, mode='full',
                                           boundary='fill')
        pyr.append(smoothed[0::2, 0::2])

    return pyr

def gen_cov(x_deriv, y_deriv, x_pos, y_pos):
    """
    Estimates the local covariance matrix of the derivative by averaging the
    outer products of the vectors x_deriv[x_pos] * x_hat and
    y_deriv[y_pos] * y_hat over a 5x5 region.

    Parameters:
    -----------

    Return:
    -------
    """
    avg = np.zeros(5, 5)
    for i in range(5):
        for j in range(5):
            xx = x_deriv[x_pos] * x_deriv[x_pos + 2 - i]
            xy = x_deriv[x_pos] * y_deriv[x_pos + 2 - j]
            yx = y_deriv[y_pos] * x_deriv[x_pos + 2 - i]
            yy = y_deriv[y_pos] * y_deriv[x_pos + 2 - j]

            avg += np.array([[xx, xy],[yx, yy]])
    avg = avg / 25.0

    return avg


def extract_features(orientation_tensor):
    """
    Performs PCA on an orientation tensorself.

    Parameters:
    -----------
    orientation :

    Returns:
    --------
    energy : sum of the two eigenvalues
    orientedness : difference of eigenvalues / sum of eigenvalues
    theta : angle of dominant eigenvector
    """
    # 2 dimensional PCA for 2D derivative
    pca = PCA(n_components = 2)
    pca.fit(orientation_tensor)

    # retrieve/compute relevant information
    eigen_values = pca.explained_variance_
    eigen_vectors = pca.components_
    energy = np.sum(eigen_values)
    orientedness = (eigen_values[0] - eigen_values[1]) / energy
    theta = np.arctan(eigen_vectors[0][1] / eigen_vectors[0][0])

    return energy, orientedness, theta
