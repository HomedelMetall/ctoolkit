from ctoolkit.global_vars.ext_libs import *
from ctoolkit.global_vars.decorators import *
from ctoolkit.tools import tools

# Subclass to manage OUTCAR file i/o
class outcar():
    # Initialize. Note that mothership variable allows for
    # POSCAR-outcar objects to interact within a VASP object.
    def __init__(self, mothership):
        self.tools = tools.tools()
        self.energy = 0.0
        self.born_charges = np.zeros([1, 3, 3], dtype=float)
        self.masses = np.zeros([1], dtype=float)
        self.volume = 0.0

        self.tol = 0.01
        self.mothership = mothership
        self.is_filled = False

    # OUTCAR i/o handler.
    # Allows to tell the library which property to extract from OUTCAR.
    # Maybe it can be automatized.
    @calculate_time
    def extract_data(self, filename, struct_Property=None, POSCAR_filename = None):
        if struct_Property is None: print("Please, specify a property to extract\n"+
                                        "from outcar from available extraction routines:\n"+
                                        "- Born charges ('BC')\n"+
                                        "- Last energy ('E')\n"+
                                        "- Last volume ('V')\n"+
				                    	"- MD_trajectory ('Traj')"
                                        "4- Energy trajectory & volumes in trajectory ('AllE')\n"+
                                        "5- Trajectory volumes ('AllV')\n"+
                                        "6- Phonons ('P')"); return -1
       
        # If user wants, poscar can be associated to the outcar structure
        if not POSCAR_filename is None:
            self.mothership.POSCAR.read_poscar(POSCAR_filename)

        outcar_file = open(filename, 'r')
        outcar_content = outcar_file.readlines()
        outcar_file.close()
        if struct_Property == 'BC':
            self.born_charges, self.masses = self.read_born_OUTCAR(outcar_content)

        elif struct_Property == 'E':
            self.energy = self.read_E_OUTCAR(outcar_content, 'last')

        elif struct_Property == 'V':
            self.volume = self.read_VOL_OUTCAR(outcar_content, 'last')
            if self.mothership.getStatus()['POSCAR']:
                if(np.abs(self.mothership.POSCAR.volume - self.volume) > self.tol):
                        print("WARNING!!! Possibly incoherent POSCAR and OUTCAR provided files?\n"+
                              " Cell volumes between the OUTCAR and POSCAR of the same VASP object\n"+
                              " have values exceeding %.8e Ang**3." % (self.tol))
                        print("POSCAR volume: %.4f OUTCAR volume: %.4f" % (self.mothership.POSCAR.volume, self.volume))

        elif struct_Property == 'P':
            self.phonons = self.read_phonons_OUTCAR(outcar_content)

        elif struct_Property == 'Traj':
            self.MDtraj, self.MDtrajdict = self.extract_MD_snapshots(outcar_content)
            # self.MDtraj = [epochs][energy, basis, positions, forces, stress_matrix]
            # self.MDtrajdict = {'energy', 'basis', 'frac_positions', 'forces', 'stress_matrix', 'stress_voigt'}
            # It is convenient here to attach POSCAR structures to this
            if not self.mothership.getStatus()['POSCAR'] is None:
                self.MDtrajStructs = self.tools.structs_from_MD(self.MDtrajdict, self.mothership.POSCAR.names, self.mothership.POSCAR.multiplicity)
                #self.MDtrajStructs = self.structs_from_MD()

        elif struct_Property == 'AllE':
            1/0 # Error

        elif struct_Property == 'AllV':
            1/0 # Error
        
        else: 
            print('Tag for physical property extraction not recognised. Check documentation (if it exists).')
            sys.exit()

        self.is_filled = True

    # Function to read BORN EFFECTIVE CHARGES from OUTCAR
    @calculate_time
    def read_born_OUTCAR(self, content, masses=None):
        Blines = content
        BDim = len(Blines)

        # MASSES WILL BE NONE UNTIL WE EXTRACT THEM FROM OUTCAR. Clean strategy.
        m_title = -1
        for i in range(1,BDim):
            if("BORN EFFECTIVE CHARGES" in Blines[i]):
                m_title = i

        m_title = m_title+2
        l = 0
        born = []
        added_rows = 0
        born_now = []
        if(self.mothership.POSCAR.is_filled == False):
            print("Please, load the associated POSCAR file of the BORN calculation to the VASP structure.")
            sys.exit()
        num_atoms = len(self.mothership.POSCAR.at_frac)
        for i in range(m_title, m_title+4*num_atoms):
            if(((i-m_title)%4)==0):
                continue
            if(((i-m_title)%4)==1):
                #if masses is None: fopen2.writelines(str(l+1) + '\t' + str(1.0) + "\n")
                #else:
                #    fopen2.writelines(str(l+1) + '\t' + str(masses[l]) + "\n")
                l = l+1
            born_now.append(np.array(Blines[i].split()[1:], dtype=float))
            added_rows += 1
            if(added_rows == 3):
                born.append(np.copy(np.transpose(np.array(born_now, dtype=float))))
                born_now = []
                added_rows = 0

        return np.copy(np.array(born, dtype=float)), np.zeros([l], dtype=float)

    # Function to read the energy from OUTCAR
    @calculate_time
    def read_E_OUTCAR(self, content, switch='last'):
        e_raw = []
        for line in content:
            if "energy  without entropy=" in line:
                e_raw.append(float(line.split()[3]))
        e_raw = np.array(e_raw, dtype=float)
    
        if switch == 'last':
            return e_raw[-1]
        if switch == 'all':
            return e_raw

    # Function to read Volumes from OUTCAR
    @calculate_time
    def read_VOL_OUTCAR(self, content, switch='last'):
        vol = []
        for i, line in enumerate(content):
            if "VOLUME and BASIS-vectors are now :" in line:
                vol.append(content[i+3].split()[4])
        vol = np.array(vol, dtype=float)

        if switch == 'last':
            return vol[-1]
        if switch == 'all':
            return vol

    # Function to read PHONONS from OUTCAR
    def read_phonons_OUTCAR(self, content):
        # Reads phonons in cartesian coordinates
        phonons = []
        end = False
        for i, line in enumerate(content):
            found = False
            if "f  =" in line or "f/i=" in line:
                found = True
                eigvals = float(line.split()[-2])
                if "f/i=" in line: eigvals = -eigvals
                phon_cart = []
                at_cart = []
                at = []
                phon = []
                for j,subline in enumerate(content[i+2:]):
                    if "--------------------------------------------------------------------------------------------------------" in content[i+j+5]:
                        end = True
                        break
                    if "f  =" in content[i+j+3] or "f/i=" in content[i+j+3]:
                        break
                    at.append(np.array(subline.split()[0:3], dtype=float))
                    phon.append(np.array(subline.split()[3:6], dtype=float))
                    #print(j, at)

            if found == True:
                phonons.append([eigvals, np.array(phon)])
            if end == True: break
        return phonons

    def extract_MD_snapshots(self, content):
        Flines = content
        FDim = len(Flines)

        # Find the positions in file of the force calculations
        for i in range(1,FDim):
            if("NIONS" in Flines[i].split()):
                numatoms = int(Flines[i].split()[11])

        i_f = 0
        for i in range(1,FDim):
            if(Flines[i]==" POSITION                                       TOTAL-FORCE (eV/Angst)\n"):
                i_f = i_f+1

        f_title = [0]*(i_f)
        i_f = 0
        for i in range(1,FDim):
            if(Flines[i]==" POSITION                                       TOTAL-FORCE (eV/Angst)\n"):
                f_title[i_f] = i
                i_f = i_f+1

        # Find the basis vectors
        s_f = 0
        for i in range(1,FDim):
            if(Flines[i]==" VOLUME and BASIS-vectors are now :\n"):
                s_f = s_f+1

        s_title = [0]*(s_f)
        s_flag = True
        s_f = 0
        for i in range(1,FDim):
            if(Flines[i]==" VOLUME and BASIS-vectors are now :\n"):
                s_title[s_f] = i
                s_f = s_f+1

        stress_f = 0
        for i in range(1,FDim):
            if(Flines[i]=="  FORCE on cell =-STRESS in cart. coord.  units (eV):\n"):
                stress_f = stress_f+1

        if (stress_f == 0):
            s_flag = False

        stress_title = [0]*(stress_f)
        stress_f = 0
        for i in range(1,FDim):
            if(Flines[i]=="  FORCE on cell =-STRESS in cart. coord.  units (eV):\n"):
                stress_title[stress_f] = i
                stress_f = stress_f+1

        counter = 0

        MDtraj = []
        MDtrajdict = {'energy':[], 'basis':[],'frac_positions':[],'forces':[],'stress_matrix':[],'stress_voigt':[]}
        for i in range(len(f_title)):
            stress = np.zeros([6])

            energy = float(Flines[f_title[i] + 2 + numatoms + 10].split()[4])
            B2 = np.zeros([3,3])
            for k in range(3):
                for j in range(3):
                    B2[j,k] = float(Flines[s_title[i] + 5 + j].split()[k])

            positions = []
            forces = []
            for k in range(numatoms):
                positions.append([float(x) for x in Flines[f_title[i] + 2 + k].split()[0:3]])
                forces.append([float(x) for x in Flines[f_title[i] + 2 + k].split()[3:6]])

            stress = [float(x) for x in Flines[stress_title[i] + 14].split()[2:2+6]]
            
            stress_voigt = np.zeros(6)
            for i in range(3):
                stress_voigt[i] = stress[i]
                stress_voigt[3+i] = stress[-i]

            stress_M = [[stress[0], stress[5], stress[4]],
            [stress[5], stress[1], stress[3]],
            [stress[4], stress[3], stress[2]]]

            MDtraj.append([[energy], B2, positions, forces, stress_M])
            MDtrajdict['energy'].append(energy)
            MDtrajdict['basis'].append(B2)
            MDtrajdict['frac_positions'].append(positions)
            MDtrajdict['forces'].append(forces)
            MDtrajdict['stress_matrix'].append(stress_M)
            MDtrajdict['stress_voigt'].append(stress_voigt)

        return MDtraj, MDtrajdict

    def structs_from_MD(self, MDtrajdict):
        num_snaps = len(MDtrajdict['energy'])
        VASPstructs = []
        for i in range(num_snaps):
            struct = self.VASP()
            # Fill poscar
            B = MDtrajdict['basis'][i]
            alat = 1.0
            names = self.mothership.POSCAR.names
            mult = self.mothership.POSCAR.multiplicity
            atoms_frac = MDtrajdict['frac_positions'][i]
            struct.POSCAR.generate_poscar(B, alat, names, mult, atoms_frac, system_name='Generated by the Crystal Toolkit')
            VASPstructs.append(struct)

        return VASPstructs
