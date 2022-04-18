# Declare exported functions
__all__ = [ 'getpath', 'get_slice']

from pathlib import Path

def getpath(*a):
    ''' Get the directory in which the current file lives 
    '''
    # Package root
    d = Path(__file__).parent.resolve()
    return d.joinpath(*a)

def get_slice(a, dim):
    ''' Returns a sliced array, given an array a. 
    dim should be a string which is either 'x', 'y' or 'z'.
    '''
    if dim=='x':
        return a[: , 0]
    elif dim=='y':
        return a[: , 1]
    elif dim=='z':
        return a[:, 2]
    else:
        # for the option dim='xyz'
        return a 