__all__ = ["tools", "decorators"]

import time
from toolkit.global_vars.timers_dict import *
from toolkit.global_vars.version import *
from toolkit.tools import * 
from toolkit.VASP import *
from toolkit.LAMMPS import *
from toolkit.GROMACS import *
from toolkit.xyz import *
from toolkit.SCALEUP import *
from toolkit.decorators import *

timers_dict['Total run'] = time.time()

from toolkit.info import *

#info()
if __name__=='__main__':
    print("HOLA")
    info()
