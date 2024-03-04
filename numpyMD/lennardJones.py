import numpy as np


def lennardJones(nat, positions):
    energy = 0.0
    forces = np.zeros(positions.shape)
    for i in range(nat):
        for j in range(i):
            dr = positions[i,:] - positions[j,:]
            dd = np.sum(dr**2)
            dd2 = 1.0 / dd
            dd6 = dd2 * dd2 * dd2
            dd12 = dd6 * dd6
            energy += 4.0 * (dd12 - dd6)
            tt = 24.0 * dd2 * (2.0 * dd12 - dd6)
            t = dr * tt
            forces[i,:] += t
            forces[j,:] -= t
    return energy, forces

