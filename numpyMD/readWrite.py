import numpy as np

def readXyz(filename,):
    #TODO: write a function that opens a file
    # 1. open the file in the read mode
    
    # 2. read the first line to get the number of atoms
    
    # 3. read the 3rd line and skip it you won't need it
    
    # 4. create a numpy array in the shape (numberOfAtoms, 3) for the positions
    
    # 5. loop over all atoms and read line by line the coordinates
    
    # 6. close the file
    
    return nat, positions


def writeXyz(filename, nat, positions, append = True):
    msg = str(nat) + '\n'
    msg += '   ' + '\n'
    for i in range(nat):
        x = positions[i,0]
        y = positions[i,1]
        z = positions[i,2]
        msg += 'LJ   ' + str(x) + '   ' + str(y) + '   ' + str(z) + '\n'

    if append:
        mode = 'a'
    else:
        mode = 'w'
    f = open(filename, mode)
    f.write(msg)
    f.close()
    return None




