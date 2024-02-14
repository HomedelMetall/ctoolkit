from toolkit.global_vars.ext_libs import *
from toolkit.global_vars.decorators import *
from toolkit.tools import tools
class GROMACS():
    def __init__(self):
        pass
    
    def __init__(self):
        self.tools = tools()
        self.type = 'GROMACS'
    
    # This will read atomic positions of an MD run. Probably easy to tune it to read velocities.
    @calculate_time
    def read_gro_trajectory(self, filename, stype=None):
        with open(filename, 'r') as f:
            # Original GROMACS units are nm, ps, C, K, kJ/mol, bar, nm/ps
            nmtoA = 10
            nmpstoAfs = 0.01
            lines = f.readlines()
            import gc
            #for i, line in enumerate(f):
            #    if(i==1):
            #        numatoms = int(line.split()[0])
            #        break
            #numsteps = int(sum([1 for _ in f])/(numatoms+3))
            numatoms = int(lines[1].split()[0])
            numsteps = int(len(lines)/(numatoms+3))
            self.atposcart = np.zeros([numsteps,numatoms,3], dtype=float)
            self.atposfrac = np.zeros([numsteps,numatoms,3], dtype=float)
            self.atvelcart = np.zeros([numsteps,numatoms,3], dtype=float)
            self.atvelfrac = np.zeros([numsteps,numatoms,3], dtype=float)
            self.index = np.zeros([numsteps, numatoms], dtype=int)
            self.atom_names = np.zeros([numsteps, numatoms], dtype=str)
            self.boxes = np.zeros([numsteps, 3, 3], dtype=float)
            self.cellA = np.zeros([numsteps], dtype=float)
            self.cellB = np.zeros([numsteps], dtype=float)
            self.cellC = np.zeros([numsteps], dtype=float)
            self.vol = np.zeros([numsteps], dtype=float)
            pattern = r'[0-9]'
            print("Steps to be read:", numsteps)
            for i in range(numsteps-1, -1, -1):
                # First, read box
                bl = np.array([float(x) for x in lines[(numatoms+3)*(i+1) - 1].split()])
                self.boxes[i] = np.transpose(np.array([[bl[0], bl[3], bl[4]],
                                [bl[6], bl[1], bl[5]],
                                [bl[7], bl[8], bl[2]]], dtype=float))*nmtoA
                
                self.cellA[i] = np.sqrt(np.dot(self.boxes[i,:,0], self.boxes[i,:,0]))
                self.cellB[i] = np.sqrt(np.dot(self.boxes[i,:,1], self.boxes[i,:,1]))
                self.cellC[i] = np.sqrt(np.dot(self.boxes[i,:,2], self.boxes[i,:,2]))
                self.vol[i] = np.linalg.det(self.boxes[i]) 

                for at in range(numatoms-1, -1, -1):
                    l = lines[2+at + (numatoms+3)*i].split()
                    if(stype=='vel'):
                        self.atvelcart[i,at,:] = np.array([float(x) for x in l[6:9]], dtype=float)*nmpstoAfs
                    if(stype=='pos'):
                        self.atposcart[i,at,:] = np.array([float(x) for x in l[3:6]], dtype=float)*nmtoA
                    if(stype==None):
                        self.atposcart[i,at,:] = np.array([float(x) for x in l[3:6]], dtype=float)*nmtoA
                        self.atvelcart[i,at,:] = np.array([float(x) for x in l[6:9]], dtype=float)*nmpstoAfs
                    self.index[i, at] = int(l[2])
                    self.atom_names[i, at] = re.sub(pattern, '', l[1]) 
                    
                # Check if we're calling well cart_to_frac
                self.atposfrac[i] = self.tools.cart_to_frac(self.boxes[i], self.atposcart[i])
                self.atvelfrac[i] = self.tools.cart_to_frac(self.boxes[i], self.atvelcart[i])
                del lines[(2+numatoms+3)*i:]
                gc.collect()
