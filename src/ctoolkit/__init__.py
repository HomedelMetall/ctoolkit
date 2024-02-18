#__all__ = ["tools", "decorators"]

import time
from ctoolkit.global_vars.timers_dict import *
from ctoolkit.global_vars.version import *
from ctoolkit.tools import * 
from ctoolkit.VASP import *
from ctoolkit.LAMMPS import *
from ctoolkit.GROMACS import *
from ctoolkit.xyz import *
from ctoolkit.SCALEUP import *
from ctoolkit.decorators import *

timers_dict['Total run'] = time.time()

from ctoolkit.info import *

#info()
if __name__=='__main__':
    pass
    #print("HOLA")
    #info()
