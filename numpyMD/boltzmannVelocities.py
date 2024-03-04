import numpy as np

def getBoltzmannVelocities(nat, positions, masses, temperature):
    # TODO: create a function that normally distributes the initial velocities
    kB = 8.617333262e-5
    # 1. draw random numbers with the standard distributions of numpy (https://numpy.org/doc/stable/reference/random/generated/numpy.random.standard_normal.html) in the shape (numberOfAtoms, 3)
    
    # 2. calculate the kinetic energy of the drawn random numbers (kineticEnergy = 0.5 * sum(v**2))
     
    # 3. calculate the target kinetic energy (targetKineticEnergy = 0.5 * kB * T)

    # 4. calculate the scaling factor for the velocities (scalingFactor = sqrt(targetKineticEnergy/kineticEnergy))
    
    # 5. scale the velocities with the scaling factor

    velocities = elim_moment(velocities)
    velocities = elim_torque(velocities, positions, masses)

    return velocities



def elim_moment(velocities):
    """
    Elimination of the momentum in the velocities
    """
    # eliminiation of momentum
    _s = np.sum(velocities, axis=0) / velocities.shape[0]
    no_moment_velocities = velocities - _s
    return no_moment_velocities

def elim_torque(velocities, positions, masses):
    """
    Elimination of the torque in the velocites
    """

    # Calculate center of mass and subtract it from positions
    total_mass = np.sum(masses)
    masses_3d = masses#np.vstack([masses] * 3).T
    weighted_positions = positions * masses_3d

    cm = np.sum(weighted_positions, axis=0) / total_mass
    weighted_positions -= cm

    # Calculate moments of inertia with both functions
    evaleria, teneria = moment_of_inertia(positions, masses)

    # New vrot calculation: Vectorized operation replacing the loop
    teneria_reshaped = teneria.T.reshape(3, 1, 3)
    vrot = np.cross(teneria_reshaped, positions[None, :, :]).transpose(1, 2, 0)

    # flatten velocities and reshape vrot to match dimensions
    velocities = velocities.flatten()
    vrot = vrot.reshape((positions.shape[0] * 3, 3), order="C")

    # normalize vrot using a loop
    for i, vec in enumerate(vrot.T):
        vrot[:, i] = normalize(vec)
        weighted_positions += cm

    # New Implementation: Vectorized operation replacing the above for loop
    # mask for elements of evaleria that are greater than 1e-10
    mask = np.abs(evaleria) > 1e-10

    # calculate alpha and update velocities using np.einsum for dot product and updating velocities
    alpha = np.einsum('ij,i->j', vrot[:, mask], velocities)
    velocities -= np.einsum('ij,j->i', vrot[:, mask], alpha)

    # reshape velocities back to the original shape
    velocities = velocities.reshape((positions.shape[0], 3))


    return velocities

def moment_of_inertia(positions, masses):
    '''
    Calcualtion of the eigenvalues and eigenvectors of the inertia tensor
    '''
    inertia_tensor = np.zeros((3, 3))
    for at, mass in zip(positions, masses):
        inertia_tensor[0, 0] += (at[1] ** 2 + at[2] ** 2)
        inertia_tensor[1, 1] += (at[0] ** 2 + at[2] ** 2)
        inertia_tensor[2, 2] += (at[0] ** 2 + at[1] ** 2)
        inertia_tensor[0, 1] -= (at[0] * at[1])
        inertia_tensor[0, 2] -= (at[0] * at[2])
        inertia_tensor[1, 2] -= (at[1] * at[2])

    inertia_tensor[1, 0] = inertia_tensor[0, 1]
    inertia_tensor[2, 0] = inertia_tensor[0, 2]
    inertia_tensor[2, 1] = inertia_tensor[1, 2]

    eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
    return eigenvalues, eigenvectors


def normalize(v):
    """
    Function that normalized a vector of arbitrary length
    """
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


