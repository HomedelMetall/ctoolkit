from toolkit.global_vars.ext_libs import *
from toolkit.global_vars.decorators import *

# Subclass to manage OUTCAR file i/o
class outcar():
    # Initialize. Note that mothership variable allows for
    # POSCAR-outcar objects to interact within a VASP object.
    def __init__(self, mothership):
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
    def extract_data(self, filename, struct_Property=None):
        if struct_Property is None: print("Please, specify a property to extract\n"+
                                        "from outcar from available extraction routines:\n"+
                                        "1- Born charges ('BC')\n"+
                                        "2- Last energy ('E')\n"+
                                        "3- Last volume ('V')\n"+
                                        "4- Energy trajectory & volumes in trajectory ('AllE')\n"+
                                        "5- Trajectory volumes ('AllV')\n"+
                                        "6- Phonons ('P')"); return -1
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
