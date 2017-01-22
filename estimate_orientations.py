import numpy as np
from project_fst import generate_rotation_matrix


def generate_noise(image):
    """Method to generate noisy images
    """
    img = image
    noise = np.zeros(size(image));
    image = img + noise
    return image

def compute_common_lines(F):
    """Method to compute the common lines 
   
    @param: A NXNX2 array of common lines where a certain l[i,j] is the common line between two planes i and j in a random image
    return: F, which a list of orientations
    """
    N = len(F)
    l = np.zeros((N, N, 2))
    for i in range(N):
        for j in range(i + 1, N):
            I_ij = np.cross(F[i][:,2], F[j][:,2])
            I_ij /= np.linalg.norm(I_ij)
            l[i,j] = np.dot((F[i][:,0:2]).T, I_ij)
            l[j,i] = np.dot((F[j][:,0:2]).T, I_ij)

    return l

def estimate_orientations(l):
    """Method to estimate orientations from the common lines generated above
        
    @param: A NXNX2 array of common lines where a certain l[i,j] is the common line between two planes i and j in a random image
    return: F, which a list of orientations
    """

    N = l.shape[0]
    for i in range(N):
        a = [generate_rotation_matrix() for i in range(10)]

    
    for rot_2, em_2 in [ (1,1), (1,-1), (-1, 1), (-1, -1)]:
        iota_2 = np.linalg.inv( np.array([v_21, em_2*w_21]).T)
        F_2[:,0:2] = np.dot( S_2, np.dot(R_2, iota_2))
        F_2[:,2] = np.cross(F_2[:,0], F_2[:,1])

        
        for rot_3, em_3 in [(1,1),(1,-1),(-1,1),(-1,-1)]:
            iota_3 = np.linalg.inv( np.array( [ v_31, em_3*w_31]).T)
            F_3[:,0:2] = np.dot( S_3, np.dot( R_3, iota_3 ))
            F_3[:,2] = np.cross(F_3[:,0], F_3[:,1])

            if np.allclose(np.dot(F_2[:,0:2], v_23), np.dot(F_3[:,0:2], v_32)):
                F[0], F[1], F[2] = F_1, F_2, F_3

    for i in range(3,N):
        V_1i = np.dot(F[0][:,0:2], l[0,i])
        V_2i = np.dot(F[1][:,0:2], l[1,i])

        v_i1 = l[i,0]
        v_i2 = l[i,1]

        iota_i = np.linalg.inv( np.array([ v_i1, v_i2 ]).T )
        R = np.array([V_1i, V_2i]).T
        F.append(np.dot(R, iota_i))

    return F


        