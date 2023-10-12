from toolkit.global_vars.version import *

class info:
    def __init__(self):
        print(self.say_hi())

    def say_hi(self):
        #os.system('clear')
        s = '\n\n\n\n'
        s += 'The Crystal Toolkit\n'
        s += 'Version %s (%s)\n' % (tk_version, tk_date)
        s += '========================================\n'
        s += 'By C. Escorihuela-Sayalero and Claudio Cazorla - carlos.escorihuela@upc.edu\n'
        s += 'Developed in the Group of Characterization of Materials of UPC\n\n'
        s += '"VASP" class:\n'
        s += ' - "POSCAR" subclass: crystal structure parameters + load functions\n'
        s += ' - "outcar" subclass: OUTCAR file process and data extraction\n\n'
        s += '"xyz" class:\n'
        s += ' - crystal structure parameters\n\n'
        s += '"tools" class:\n'
        s += ' - Cell operations (from structures or structure classes, interpolation...) \n'
        s += ' - Structure class management operations (copy structures, transform types...)\n'
        s += ' - 3rd Order Birch-Murnaghan Equation of State (from VASP calculations)\n'
        s += ' - Implementation of E(V, params) and P(V, params) equations\n'
        s += ' - V(P) numerical solution through bisect method\n'
        s += ' - Automation: mass-load of VASP structures and OUTCAR files\n'
        s += ' - Random tools: file I/O handlers, utilities for matrix manipulations, etc\n\n'

        s += '\n\n Below, your script run:\n'
        s += ' -----------------------------------\n'

        return s
