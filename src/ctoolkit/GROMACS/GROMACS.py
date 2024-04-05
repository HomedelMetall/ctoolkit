from ctoolkit.global_vars.ext_libs import *
from ctoolkit.global_vars.decorators import *
from ctoolkit.tools import tools
import pickle
class GROMACS():
    def __init__(self):
        self.tools = tools.tools()
        self.type = 'GROMACS'
    
    # This will read atomic positions of an MD run. Probably easy to tune it to read velocities.
    @calculate_time
    def read_gro_trajectory(self, filename, stype=None):
        # Original GROMACS units are nm, ps, C, K, kJ/mol, bar, nm/ps
        nmtoA = 10
        nmpstoAfs = 0.01
        with open(filename, 'r') as f:
            ct = 0
            for l in f:                
                ct += 1

            numlines = ct
            f.seek(0,0)
            f.readline()
            numatoms = int(f.readline().split()[0])
            numsteps = int(numlines/(numatoms+3))
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
            f.seek(0, 0) 
            for i in range(numsteps):
                f.readline()
                f.readline()
                for at in range(numatoms):
                    l = f.readline().split()
                    self.atposcart[i,at,:] = np.array([float(x) for x in l[3:6]], dtype=float)*nmtoA
                    self.atvelcart[i,at,:] = np.array([float(x) for x in l[6:9]], dtype=float)*nmpstoAfs
                    self.index[i, at] = int(l[2])
                    self.atom_names[i, at] = re.sub(pattern, '', l[1])

                bl = np.array([float(x) for x in f.readline().split()])
                self.boxes[i] = np.transpose(np.array([[bl[0], bl[3], bl[4]],
                                [bl[6], bl[1], bl[5]],
                                [bl[7], bl[8], bl[2]]], dtype=float))*nmtoA

                self.cellA[i] = np.sqrt(np.dot(self.boxes[i,:,0], self.boxes[i,:,0]))
                self.cellB[i] = np.sqrt(np.dot(self.boxes[i,:,1], self.boxes[i,:,1]))
                self.cellC[i] = np.sqrt(np.dot(self.boxes[i,:,2], self.boxes[i,:,2]))
                self.vol[i] = np.linalg.det(self.boxes[i])

                # Check if we're calling well cart_to_frac
                self.atposfrac[i] = self.tools.cart_to_frac(self.boxes[i], self.atposcart[i])
                self.atvelfrac[i] = self.tools.cart_to_frac(self.boxes[i], self.atvelcart[i])
                sys.stdout.write("%d/%d\r" % (i, numsteps))
                sys.stdout.flush()

            """
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
                    ''' A slow selector
                    if(stype=='vel'):
                        self.atvelcart[i,at,:] = np.array([float(x) for x in l[6:9]], dtype=float)*nmpstoAfs
                    if(stype=='pos'):
                        self.atposcart[i,at,:] = np.array([float(x) for x in l[3:6]], dtype=float)*nmtoA
                    if(stype==None):
                        self.atposcart[i,at,:] = np.array([float(x) for x in l[3:6]], dtype=float)*nmtoA
                        self.atvelcart[i,at,:] = np.array([float(x) for x in l[6:9]], dtype=float)*nmpstoAfs
                    '''
                    self.atposcart[i,at,:] = np.array([float(x) for x in l[3:6]], dtype=float)*nmtoA
                    self.atvelcart[i,at,:] = np.array([float(x) for x in l[6:9]], dtype=float)*nmpstoAfs
                    self.index[i, at] = int(l[2])
                    self.atom_names[i, at] = re.sub(pattern, '', l[1]) 
                    
                # Check if we're calling well cart_to_frac
                self.atposfrac[i] = self.tools.cart_to_frac(self.boxes[i], self.atposcart[i])
                self.atvelfrac[i] = self.tools.cart_to_frac(self.boxes[i], self.atvelcart[i])
                del lines[(2+numatoms+3)*i:]
                sys.stdout.write("%d/%d\r" % (numsteps-i, numsteps))
                sys.stdout.flush()
                #gc.collect()

            """

    def __getstate__(self):
        state = self.__dict__.copy()
        for i, st in enumerate(state):
            if not self.is_picklable(st):
                print("Not picklable", i, state)
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def test_picklability(self, obj, name="Object"):
        try:
            pickle.dumps(obj)
            print(f"{name} is picklable.")
        except (pickle.PicklingError, TypeError) as e:
            print(f"{name} is NOT picklable. Reason: {e}")

    def check_instance_picklability(self, instance):
        for attr_name in dir(instance):
            # Filter out special and private attributes
            if not attr_name.startswith("__"):
                attr = getattr(instance, attr_name)
                self.test_picklability(attr, name=attr_name)
    
    def is_picklable(self, obj):
        try:
            pickle.dumps(obj)
            return True
        except (pickle.PicklingError, TypeError):
            return False

    @calculate_time
    def read_gro_trajectory_parallel(self, filename):
        with open(filename, 'r') as f:
            # Original GROMACS units are nm, ps, C, K, kJ/mol, bar, nm/ps
            nmtoA = 10
            nmpstoAfs = 0.01
            lines = f.readlines()
            numatoms = int(lines[1].split()[0])
            numsteps = int(len(lines)/(numatoms+3))
            
            print("Steps to be read:", numsteps)

        with self.tools.parallel_manager() as m:
            boxes = m.list([0]*numsteps)
            cellA, cellB, cellC = m.list([.0]*numsteps), m.list([.0]*numsteps), m.list([.0]*numsteps)
            vol = m.list([.0]*numsteps)
            atom_names = m.list([''*numatoms]*numsteps)
            index = m.list([0]*numsteps)
            atposfrac = m.list([0]*numsteps)
            atvelfrac = m.list([0]*numsteps)
            atposcart = m.list([0]*numsteps)
            atvelcart = m.list([0]*numsteps)
            
            # PARALLEL RUN
            global lock
            lock = self.tools.get_lock()
            pool = self.tools.create_pool(lock)
            targs = [(step_id, numatoms, lines, boxes, cellA, cellB, cellC,
                    vol, atom_names, index,
                    atposcart, atposfrac,
                    atvelcart, atvelfrac, ) for step_id in range(numsteps)]
           
            res = self.tools.run_parallel(read_step, targs, chksz=30, pool=pool)
            self.tools.close_pool(pool)
            

            self.boxes = np.array(boxes)
            self.cellA, self.cellB, self.cellC = np.array(cellA), np.array(cellB), np.array(cellC)
            self.vol = np.array(vol)
            self.atom_names = np.array(atom_names, dtype=str)
            self.index = np.array(index, dtype=int)
            self.atposcart = np.array(atposcart)
            self.atvelcart = np.array(atvelcart)

            self.atposfrac = np.copy(self.atposcart)
            self.atvelfrac = np.copy(self.atvelcart)
           
            for i in range(len(self.atposcart)):
                self.atposfrac[i] = self.tools.cart_to_frac(self.boxes[i], self.atposcart[i])
                self.atvelfrac[i] = self.tools.cart_to_frac(self.boxes[i], self.atvelcart[i])

            del lines

## PARALLEL SECTION ##
# When you're coding your thing and calling ctoolkit
# to instantiate a parallel job, all is smooth.
# However, from inside, it is very difficult
# to establish the notion of SELF.
# So if we try to create looping functions
# using self. methods, multiprocessing gets
# crazy because can't pickle the objects.

# INSTEAD.

# We do it here.
# We do it nasty.
# We do it *OUTSIDE*
def read_step(step_id, numatoms, lines,
                    _boxes, _cellA, _cellB, _cellC,
                    _vol, _atom_names, _index,
                    _atposcart, _atposfrac,
                    _atvelcart, _atvelfrac) -> None:
    i = int(step_id+0)
    nmtoA = 10
    nmpstoAfs = 0.01
    bl = np.array([float(x) for x in lines[(numatoms+3)*(i+1) - 1].split()])

    boxes = np.transpose(np.array([[bl[0], bl[3], bl[4]],
                           [bl[6], bl[1], bl[5]],
                           [bl[7], bl[8], bl[2]]], dtype=float))*nmtoA

    cellA = np.sqrt(np.dot(boxes[:,0], boxes[:,0]))
    cellB = np.sqrt(np.dot(boxes[:,1], boxes[:,1]))
    cellC = np.sqrt(np.dot(boxes[:,2], boxes[:,2]))
    vol = np.linalg.det(boxes)

    atposcart = np.zeros([numatoms, 3], dtype=float)
    atvelcart = np.zeros([numatoms, 3], dtype=float)
    index = np.zeros([numatoms], dtype=float)
    atom_names = ['']*numatoms#np.zeros([numatoms], dtype=str)
    pattern = r'[0-9]'
    for at in range(numatoms):
        l = lines[2+at + (numatoms+3)*i].split()
        atposcart[at,:] = np.array([float(x) for x in l[3:6]], dtype=float)*nmtoA
        atvelcart[at,:] = np.array([float(x) for x in l[6:9]], dtype=float)*nmpstoAfs
        index[at] = int(l[2])
        atom_names[at] = re.sub(r'\d.*', '', l[1])

    # Split writeout to memory so as to use a datalock
    with lock:
        _boxes[step_id] = boxes
        _cellA[step_id], _cellB[step_id], _cellC[step_id] = cellA, cellB, cellC
        _vol[step_id] = vol
        _atom_names[step_id] = atom_names
        _index[step_id] = index
        _atposcart[step_id] = np.array(atposcart)
        _atvelcart[step_id] = np.array(atvelcart)

