#import statements for various packages needed to compute values and plot end result
import numpy as np
import matplotlib as plt
from scipy.interpolate import RegularGridInterpolator
from MRCFile import MRCFile
import matplotlib.pyplot as plt


def n_matrix(n, range):
    """
    @param n -  molecule's shape
    @param range: range value to create a N sized matrix to be used for the image plane space
    return: a numpy ndarray that represents the image projected using mol and R
    """
    return np.concatenate((range.reshape(n, 1, 1, 1) * np.ones(n * n).reshape(1, n, n, 1), range.reshape(1, n, 1, 1) * np.ones(n * n).reshape(n, 1, n, 1),
                          range.reshape(1, 1, n, 1) * np.ones(n * n).reshape(n, n, 1, 1)), axis=3)



def project_fst(mol, R):
    """
    @param mol -  We pass in an NxNxN array that contains the values for rho
    @param R: a rotation matrix
    return: a numpy ndarray that represents the image projected using mol and R
    """
    N = mol.shape[0]
    interpolation = RegularGridInterpolator((np.linspace(-1, 1, N), np.linspace(-1, 1, N), np.linspace(-1, 1, N)), mol, method='linear', bounds_error=False, fill_value=0)
    image = np.real(np.fft.ifft2(np.fft.fftn(interpolation(np.dot(n_matrix(N, np.linspace(-1, 1, N)), R)))[:, :, 0]))
    # n_matrix returns a matrix of size N that returns a value that can be transformed into the image plane space
    return image

def generate_rotation_matrix():
    """
    Method that is able to randomly generate rotation matrices that we can use to test project_fst
    """
    a = 2 * np.random.rand(3) - 1 
    b = np.cross(a, 2 * np.random.rand(3) - 1) 
    c = np.cross(a, b) 
    a /= np.linalg.norm(a)
    b /= np.linalg.norm(b)
    c /= np.linalg.norm(c)
    return np.array([a, b, c]).T

molecule = MRCFile('zika_153.mrc')
#image = project_fst(molecule.data, np.array([[np.sqrt(2)/2, np.sqrt(2)/2, 0], [0, 0, 1], [np.sqrt(2)/2, -np.sqrt(2)/2, 0]]))
image = project_fst(molecule.data, generate_rotation_matrix()) #running the project_fst function
plt.imshow(image) #plotting the image to a matplotlib.pyplot that we can visualize
plt.show()
