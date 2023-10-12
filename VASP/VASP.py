from toolkit.global_vars.ext_libs import *
from toolkit.VASP.POSCAR import *
from toolkit.VASP.outcar import *

# Main structure class VASP
# Subclasses POSCAAR and outcar
# handle file i/o in crystal toolkit format.
class VASP():
    def __init__(self):
        self.POSCAR = POSCAR(self)
        self.outcar = outcar(self)
        self.type = 'VASP'

    @calculate_time
    def getStatus(self):
        return {'POSCAR' : self.POSCAR.is_filled, 'OUTCAR' : self.outcar.is_filled}

    # Override to POSCAR.write
    @calculate_time
    def write(self, filename=None):
        self.POSCAR.write(filename)
