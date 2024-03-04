import numpy as np
from ase.io import read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.neighborlist import NeighborList
from ase.calculators.calculator import Calculator, all_changes
from ase.stress import full_3x3_to_voigt_6_stress
from numba import jit, cuda



class LennardJones(Calculator):

    implemented_properties = ['energy', 'forces']
    # implemented_properties += ['stress', 'stresses']  # bulk properties
    default_parameters = {
        'epsilon': 1.0,
        'sigma': 1.0,
        'rc': None,
        'ro': None,
        'smooth': False,
    }
    nolabel = True

    def __init__(self, **kwargs):


        Calculator.__init__(self, **kwargs)

        if self.parameters.rc is None:
            self.parameters.rc = 3 * self.parameters.sigma

        if self.parameters.ro is None:
            self.parameters.ro = 0.66 * self.parameters.rc

        self.nl = None

    def calculate(
        self,
        atoms=None,
        properties=None,
        system_changes=all_changes,
    ):
        if properties is None:
            properties = self.implemented_properties

        Calculator.calculate(self, atoms, properties, system_changes)

        e, f = LennardJones.force(atoms.get_positions())

        self.results['energy'] = e
        self.results['forces'] = f

    @staticmethod
    @jit(nopython=True, parallel=False)
    def force(ats):
        nat = ats.shape[0]
        e = 0.0
        f = np.zeros(ats.shape)

        for i in range(nat):
            for j in range(i):
                dr = ats[i,:] - ats[j,:]
                dd = np.sum(dr**2)
                dd2 = 1.0 / dd
                dd6 = dd2 * dd2 * dd2
                dd12 = dd6 * dd6
                e += 4.0 * (dd12 - dd6)
                tt = 24.0 * dd2 * (2.0 * dd12 - dd6)
                t = dr * tt
                f[i,:] += t
                f[j,:] -= t
        return e, f


def main():
    # specify some parameter
    numberOfSteps = 20000
    temperature = 10000.0
    dt = 0.01
    outputFilename = 'output_structures.xyz'

    # TODO: Read xyz file (https://wiki.fysik.dtu.dk/ase/ase/io/io.html)
    inputFilename = 'lj38.xyz'
       
    # TODO: Create and attach calculator to the atoms object (https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html)
     
    # TODO: initialize velocities according to temperature (https://wiki.fysik.dtu.dk/ase/ase/md.html)
    
    # eliminate the torque and the momentum in the velocities
    velocities = structure.get_velocities()
    positions =  structure.get_positions()
    masses = structure.get_masses()
    velocities = elim_moment(velocities)
    velocities = elim_torque(velocities,positions,masses)
    structure.set_velocities(velocities)
    
    # TODO: calculate initial forces and energy (https://wiki.fysik.dtu.dk/ase/ase/atoms.html)
    
    
    log(structure, potentialEnergy)
    write(outputFilename, structure)
    # start md loop
    for step in range(numberOfSteps):
        # TODO: get positions, velocity and masses (https://wiki.fysik.dtu.dk/ase/ase/atoms.html)
        
        
        # TODO: update the positions according to the verlet algorithm (x <- x + dt*v + 0.5*dt**2*f) (https://wiki.fysik.dtu.dk/ase/ase/atoms.html)
        
        

        # TODO: update the forces (https://wiki.fysik.dtu.dk/ase/ase/atoms.html)
        

        # TODO: update the velocities (v <- v+0.5*dt*((f+f_new))) (https://wiki.fysik.dtu.dk/ase/ase/atoms.html)
        
        
        
        log(structure, potentialEnergy,step)
        if step%100 == 0:
            write(outputFilename, structure, append = True)
        
        # TODO: update the new forces to the forces (f_old <- f_new)
        
        



def log(structure, potentialEnergy,step = None):
    if step is None:
        msg = 'MD step    epot         ekin          etot'
        step = 0
        print(msg)
    kineticEnergy = getKineticEnergy(structure)     
    totalEnergy = potentialEnergy + kineticEnergy
    msg = "{:4d}      {:1.5f}      {:1.5f}       {:1.5f}".format(
        step + 1,
        potentialEnergy,
        kineticEnergy,
        totalEnergy
    )
    print(msg)
    return None


def getKineticEnergy(structure):
    velocities = structure.get_velocities()
    kineticEnergy = 0.5 * np.sum(velocities * velocities)
    return kineticEnergy



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
    weighted_positions = positions 

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



if __name__ == '__main__':
    main()
