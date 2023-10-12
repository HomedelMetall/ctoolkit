from toolkit.global_vars.ext_libs import *

# Class to manage xyz formats.
class xyz():
    # Initialization
    def __init__(self):
        self.tools = tools()
        self.title = '\n'
        self.B = np.zeros([3,3], dtype=float)
        self.a_lat = 0
        self.names = []
        self.multiplicity = [0]
        self.namelist = []
        self.atom_id = []

        self.at_frac = np.array(np.array([.0, .0, .0], dtype=float))
        self.at_cart = np.array(np.array([.0, .0, .0], dtype=float))
        self.type='xyz'

    # Function to write out the .xyz file.
    @calculate_time
    def write(self, filename=None):
        s = ''
        s += '%d\n' % (len(self.at_cart))
        s += (self.title)

        print(self.namelist)
        for atom in range(len(self.at_cart)):
            s += '%s\t%.8f\t%.8f\t%.8f' % (self.namelist[atom], self.at_cart[atom, 0], self.at_cart[atom, 1], self.at_cart[atom, 2])
            if atom < len(self.at_cart)-1 : s += '\n'

        if filename != None:
            fopen = open(filename,'w')
            fopen.write(s)
            fopen.close()
        else:
            pass
            #print(s)
        return s

    # Function to resort the atoms given a reordering list
    @calculate_time
    def resort_atoms(self, order):
        new_cart = np.zeros([len(self.at_cart), 3], dtype=float)
        new_frac = np.zeros([len(self.at_frac), 3], dtype=float)
        new_names = np.copy(self.namelist)
        
        for i, j in enumerate(order):
            new_cart[i] = self.at_cart[j]
            new_frac[i] = self.at_frac[j]
            new_names[i] = self.namelist[j]
        
        self.at_cart = np.copy(new_cart)
        self.at_frac = np.copy(new_frac)
        self.namelist = np.copy(new_names)
