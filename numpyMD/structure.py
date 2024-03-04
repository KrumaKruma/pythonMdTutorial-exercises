import numpy as np
from boltzmannVelocities import getBoltzmannVelocities
from lennardJones import lennardJones
from readWrite import writeXyz


class Structure():
    def __init__(self, nat, positions):
        self.nat = nat
        self.positions = positions
        self.masses = np.ones((nat,3))

    def initializeVelocities(self, temperature):
        velocities = getBoltzmannVelocities(self.nat, self.positions, self.masses, temperature)
        self.velocities = velocities
        return None

    def getEnergiesAndForces(self,):
        energy, forces = lennardJones(self.nat, self.positions)
        return energy, forces





