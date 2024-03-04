import numpy as np
from readWrite import readXyz, writeXyz
from structure import Structure


def main():
    # specify some parameter
    numberOfSteps = 20000
    temperature = 10000.0
    dt = 0.01
    outputFilename = 'output_structures.xyz'

    # TODO: Read xyz file
    inputFilename = 'lj38.xyz'
    
    # TODO: create structure object

    # TODO: initialize velocities according to temperature (use the structure object)

    # TODO: calculate initial forces and energy (use the structure object)
    log(structure, potentialEnergy)
    writeMD(outputFilename, structure, append = False)
    # start md loop
    for step in range(numberOfSteps):
        # TODO: get positions, velocity and masses (use the structure object)




        # TODO: update the positions (x <- x+dt*v+0.5*dt**2*f/m) and change them in the structure object



        # TODO: update the forces (use the structure object)
        

        # TODO: update the velocities (v <- v+0.5*dt((f+f_new)/masses)) and change them in the structure object
        
        log(structure, potentialEnergy,step)
        if step%100 == 0:
            writeMD(outputFilename, structure, append=True)
        
        # TODO: update the forces (f <- f_new)


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
    masses = structure.masses
    velocities = structure.velocities
    kineticEnergy = 0.5 * np.sum(masses * velocities * velocities)
    return kineticEnergy


def writeMD(outputFilename, structure, append):
    nat = structure.nat
    positions = structure.positions
    writeXyz(outputFilename,nat,positions,append=append)
    return None



if __name__ == '__main__':
    main()
