# Declare exported functions
__all__ = [ 'MSDparams', 'VACFparams']

from dataclasses import dataclass, field

@dataclass
class MSDparams:
    ''' Holds various parameters required to perform an MSD calculation.
    These can be read in from a TOML file  
    '''
    trajectory: str
    trajectory_file_type: str = None 
    max_lag_time: int = None 
    first_time_origin: int = 0
    first_lag_time: int = 1
    step_size_lag_time: int = 1  
    dimension: str = 'all'
    use_center_of_mass: bool = False 


@dataclass
class VACFparams:
    ''' Parameters required for a VACF calculation. 
    Can be read in from a TOML file. 
    Currently, the velocity averaged over the atoms in the molecules
    is hard-coded for water (TODO: make general). 

    velocity_center_of_mass: If this is True, then the velocity 
    of the COM of the molecule is calculated
    '''
    trajectory: str
    trajectory_file_type: str = None 
    max_lag_time: int = None 
    first_time_origin: int = 0
    first_lag_time: int = 1
    step_size_lag_time: int = 1  
    velocity_center_of_mass: bool = False 