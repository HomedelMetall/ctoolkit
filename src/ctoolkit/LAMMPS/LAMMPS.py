from ctoolkit.global_vars.ext_libs import *
from ctoolkit.global_vars.decorators import *
from ctoolkit.tools import tools

# Class to manage LAMMPS file i/o
class LAMMPS():
    # Initialijze
    def __init__(self):
        self.tools = tools.tools()
        pass

    # This will read atomic positions of an MD run. Probably easy to tune it to read velocities.
    @calculate_time
    def read_snapshots_func(self, filename, atomsline):
        fopen = open(filename, 'r')
        lines = fopen.readlines()
        fopen.close()

        numatoms = int(lines[3].split()[0])
        boundbox_mrk = []
        atoms_mrk = []
        nsteps = 0
        for i, line in enumerate(lines):
            if "BOX BOUNDS" in line:
                boundbox_mrk.append(i+1)
            if atomsline in line:
                atoms_mrk.append(i+1)
            if "TIMESTEP" in line:
                nsteps += 1

        # Boundbox processing
        box = np.zeros([nsteps, 3,3], dtype=float)
        for i, mrk in enumerate(boundbox_mrk):
            box[i, 0, 0] = float(lines[mrk].split()[1])-float(lines[mrk].split()[0])
            box[i, 1, 1] = float(lines[mrk+1].split()[1])-float(lines[mrk+1].split()[0])
            box[i, 2, 2] = float(lines[mrk+2].split()[1])-float(lines[mrk+2].split()[0])

        atoms = np.zeros([nsteps, numatoms, 3], dtype=float)
        id_lammps = np.zeros([nsteps, numatoms], dtype=int)
        for i, mrk in enumerate(atoms_mrk):
            for j in range(numatoms):
                if atomsline == 'ATOMS type x y z':
                    id_lammps[i,j] = int(lines[mrk+j].split()[0])
                    atoms[i,j, :] = np.array(lines[mrk+j].split()[1:4], dtype=float)
                else:
                    atoms[i,j, :] = np.array(lines[mrk+j].split()[0:3], dtype=float)

        self.id_lammps, self.boxes, self.atoms = np.copy(id_lammps), np.copy(box), np.copy(atoms)

    #This is just a wrapper
    @calculate_time
    def read_snapshots(self, filename, style='pos'):
        if style == 'typepos':
            atomsline = 'ATOMS type x y z'
        if style == 'pos':
            atomsline = "ATOMS x y z"
        if style == 'vel':
            atomsline = "ATOMS vx vy vz"

        self.read_snapshots_func(filename, atomsline)

    # A function to read LAMMPS output.
    # This is trivial to modify but for now the 
    # reader needs that the main output has the following structure:
    #
    #   "Step Temp Press Cella Cellb Cellc Volume PotEng"
    #
    @calculate_time
    def read_LAMMPS_output(self, filename):
        fopen = open(filename, 'r')
        lines = fopen.readlines()
        fopen.close()

        # We need this because LAMMPS closes/engages a new loop every time
        # a run command is issued.
        line_starts = []
        line_ends = []
        for i, line in enumerate(lines):
            if "Step" in line:            
                sample = line
                line_starts.append(i+1)
            if "Loop time" in line:
                line_ends.append(i)

        if line_ends == []: # simulation didn't finish
            line_ends.append(len(lines))

        num_loops = len(line_starts)
        # We make a generic reader using an automated dictionary
        dict_output = {}
        for element in sample.split():
            dict_output[element] = []

        data = []
        for iloop in range(num_loops):
            # Skip the first result as it's the same as the last loop last step!
            for i in range(line_starts[iloop]+1, line_ends[iloop]):
                d = []
                for j, val in enumerate(lines[i].split()):
                    dict_output[sample.split()[j]].append(float(val))

        # Convert to numpy arrays the contents of the dictionary
        for key in list(dict_output.keys()):
            dict_output[key] = np.array(dict_output[key])
        
        # The keys here are thus the ones provided by LAMMPS itself

        return dict_output
