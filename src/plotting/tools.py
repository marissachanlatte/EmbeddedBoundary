import numpy as np

def cell_average(array, n):
    ''' Turns 1d flattened array (representing a 2d mxm array)
        into an a 1d flattened array representing a 2d nxn array
        by averaging cell values.'''
    # Check m divisible by n
    m = int(np.sqrt(len(array)))
    if m%n != 0:
        raise ValueError('Array not divisible by n.')

    block_size = int(m/n)
    averaged_array = np.zeros(n**2)

    count = 0
    for j in range(n):
        for i in range(n):
            running_sum = 0
            for k in range(block_size):
                for l in range(block_size):
                    running_sum += array[block_size*i + m*block_size*j+ k + m*l]
            averaged_array[count] = running_sum/block_size**2
            count += 1
    
    return averaged_array
