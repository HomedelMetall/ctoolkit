from toolkit.global_vars.ext_libs import *
from toolkit.tools import tools
class SCALEUP():
    def __init__(self):
        #self.read_poscar()
        self.tools = tools()
        self.title = 'Toolkit_generated\n'
        self.B = np.zeros([3,3], dtype=float)
        self.SC = np.array([0,0,0], dtype=int)
        self.a_lat = 1.0
        self.nat = 0
        self.names = []
        self.namelist = []
        self.nspecies = 0
        self.atom_id = []
        self.at_type = []
        self.cell_id = []
        self.lattype = '\n'
        self.multiplicity = [0]
        self.volume = 0.0
        self.at_frac = np.array(np.array([.0, .0, .0], dtype=float))
        self.at_cart = np.array(np.array([.0, .0, .0], dtype=float))
        self.type = "SCALEUP"

    def apply_distortion(self):
        self.B = np.dot((np.eye(3)+self.strains), self.B)
        for i in range(len(self.at_frac)):
            self.at_cart[i] = np.dot(self.B, self.at_frac[i]) + self.disp[i]
            self.at_frac[i] = self.tools.cart_to_frac(self.B, self.at_cart[i])

    def read_restart_file(self, file):
        lines = open(file, 'r').readlines()
        SF = np.array([float(lines[3].split()[x]) for x in range(6)], dtype=float)
        self.strains = np.array([[SF[0], SF[5],SF[4]],
                                [SF[5], SF[1],SF[3]],
                                [SF[4], SF[3],SF[2]]], dtype=float)
        self.disp = np.zeros([self.nat, 3], dtype=float)
        for i in range(self.nat):
            self.disp[i] = np.array([float(lines[4+i].split()[5+x]) for x in range(3)])

        self.disp = 0.52917721*self.disp

    def read_structure_file(self, file):
        lines = open(file, 'r').readlines()
        self.SC = np.array([float(lines[0].split()[x]) for x in range(3)], dtype=int)
        self.nat = int(lines[1].split()[0])*np.prod(self.SC)
        self.nspecies = int(lines[1].split()[1])
        self.names = lines[2].split()[:]
        self.B = 0.52917721*np.array([[float(lines[3].split()[3*x+y]) for y in range(3)] for x in range(3)], dtype=float)
        self.at_frac = np.zeros([self.nat, 3], dtype=float)
        self.at_cart = np.zeros([self.nat, 3], dtype=float)
        self.cell_id = np.zeros([self.nat, 3], dtype=int)
        self.atom_id = np.zeros([self.nat], dtype=int)
        self.at_type = np.zeros([self.nat], dtype=int)
        self.multiplicity = np.zeros([self.nspecies], dtype=int)
        for i in range(self.nat):
            self.at_cart[i] = 0.52917721*np.array([float(lines[4+i].split()[5+x]) for x in range(3)])
            self.at_frac[i] = self.tools.cart_to_frac(self.B, self.at_cart[i])
            self.cell_id[i] = np.array([int(lines[4+i].split()[0+x]) for x in range(3)], dtype=int)
            self.atom_id[i] = int(lines[4+i].split()[3])-1
            self.at_type[i] = int(lines[4+i].split()[4])-1
            self.multiplicity[self.at_type[i]] += 1
        self.at_frac = self.at_frac
