# Declare exported functions
__all__ = [ 'MSDparams']

from dataclasses import dataclass, field

@dataclass
class MSDparams:
    trajectory: str
    trajectory_file_type: str = None 
    max_lag_time: int = None 
    first_time_origin: int = 0
    first_lag_time: int = 1
    step_size_lag_time: int = 1  
    dimension: str = 'xyz'