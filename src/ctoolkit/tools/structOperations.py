from ctoolkit.global_vars.ext_libs import *
from ctoolkit.global_vars.tk_lib_VASP import *
from ctoolkit.global_vars.tk_lib_SCALEUP import *
from ctoolkit.global_vars.tk_lib_XYZ import *
from ctoolkit.global_vars.tk_lib_GROMACS import *
class structOperations:
    @calculate_time
    def copy_structure_vasp(self, structure):
        new_VASP = VASP()
        new_POSCAR = new_VASP.POSCAR

        new_POSCAR.title = structure.POSCAR.title + ''
        new_POSCAR.B = np.copy(structure.POSCAR.B)
        new_POSCAR.a_lat = structure.POSCAR.a_lat + 0.0
        new_POSCAR.names = structure.POSCAR.names
        new_POSCAR.lattype = structure.POSCAR.lattype + ''
        new_POSCAR.multiplicity = np.copy(structure.POSCAR.multiplicity)
        new_POSCAR.at_frac = np.copy(structure.POSCAR.at_frac)
        new_POSCAR.at_cart = np.copy(structure.POSCAR.at_cart)
        new_POSCAR.atom_id = np.copy(structure.POSCAR.atom_id)
        new_POSCAR.namelist = structure.POSCAR.namelist
        new_POSCAR.type = "VASP"
        
        new_OUTCAR = new_VASP.outcar
        new_OUTCAR.energy = structure.outcar.energy + 0.0
        new_OUTCAR.born_charges = structure.outcar.born_charges + 0.0 
        new_OUTCAR.masses = structure.outcar.masses + 0.0 
        new_OUTCAR.volume = structure.outcar.volume + 0.0 

        return new_VASP

    @calculate_time
    def copy_structure_scaleup(self, structure):
        new_SCUP = SCALEUP()
        new_SCUP.title = structure.title
        new_SCUP.B = structure.B
        new_SCUP.SC = structure.SC
        new_SCUP.alat = structure.a_lat
        new_SCUP.nat = structure.nat
        new_SCUP.names = structure.names
        new_SCUP.namelist = structure.namelist
        new_SCUP.nspecies = structure.nspecies
        new_SCUP.atom_id = structure.atom_id
        new_SCUP.at_type = structure.at_type
        new_SCUP.cell_id = structure.cell_id
        new_SCUP.lattype = structure.lattype
        new_SCUP.multiplicity = structure.multiplicity
        new_SCUP.volume = structure.volume
        new_SCUP.at_frac = structure.at_frac
        new_SCUP.at_cart = structure.at_cart
        new_SCUP.type = SCALEUP

        return new_SCUP

    @calculate_time
    def copy_structure(self, structure):
        if structure.type == "VASP":
            new_structure = self.copy_structure_vasp(structure)

        if structure.type == "SCALEUP":
            new_structure = self.copy_structure_scaleup(structure)

        return new_structure 

    def lammps_to_xdatcar(self, LAMMPS_structure, output_xdat, namelist=None):
        #header = '\n'.join(poscar_str.split('\n')[:7]) # TO BE DESCRIBED FROM LAMMPS STRUCTURE, easy-peasy
        # Example: only for MAPI for now, but here is the juice
        # We make a list containing the number of different specie id's from lammps, with the same order
        # That way we can automatically count atoms
        #namelist = ['Pb', 'I', 'H', 'N', 'C', 'H3']
        #namelist = ['C', 'B']

        # Get number of different types
        typelist = []
        for _id in LAMMPS_structure.id_lammps[0]:
            if _id in typelist: continue
            typelist.append(_id)
        num_species = len(typelist)
        if namelist == None:
            namelist = []
            for i in range(num_species):
                namelist.append('A%d' % (i+1))
        else:
            if num_species != len(namelist):
                print("The number of species inputted is different than the length of the namelist!")
                sys.exit()

        reorder_namelist = range(num_species)

        sorted_namelist = []
        for i in reorder_namelist:
            sorted_namelist.append(namelist[i])
        print("Sorting XDATCAR atomic species as:")
        print(" ".join(namelist) + " -> " + " ".join(sorted_namelist))

        sorted_typelist = []
        for i in range(len(typelist)):
            sorted_typelist.append(typelist[reorder_namelist[i]])

        print("Number of species identified: %d" % (num_species))
        # Get number of atoms of each specie
        num_atoms = np.zeros([num_species], dtype=int)
        for i, _id in enumerate(typelist):
            for _idlammps in LAMMPS_structure.id_lammps[0]:
                if _idlammps == _id: num_atoms[i] += 1
        header = ''.join(sorted_namelist) + '\n'
        header += '1.0\n'
        for i in range(3):
            header += '%.8e\t%.8e\t%.8e\n' % (LAMMPS_structure.boxes[0, i, 0], LAMMPS_structure.boxes[0, i, 1], LAMMPS_structure.boxes[0, i, 2])
        header += ' '.join(sorted_namelist) + '\n'
        sorted_num_atoms = np.zeros([num_species], dtype=int)
        for i in range(len(num_atoms)):
            # Here we reorder on site the atoms count on-the-fly according to the defined order for the namelist
            header += '%d ' % (num_atoms[reorder_namelist[i]])
            sorted_num_atoms[i] = num_atoms[reorder_namelist[i]]

        print(header)
        # Print new XDATCAR
        atoms_frac = np.zeros([len(LAMMPS_structure.atoms), len(LAMMPS_structure.atoms[0]), 3])
        atoms_cart = np.zeros([len(LAMMPS_structure.atoms), len(LAMMPS_structure.atoms[0]), 3])
        sorted_atoms_frac = np.zeros([len(LAMMPS_structure.atoms), len(LAMMPS_structure.atoms[0]), 3])
        sorted_atoms_cart = np.zeros([len(LAMMPS_structure.atoms), len(LAMMPS_structure.atoms[0]), 3])
        for i in range(len(LAMMPS_structure.atoms)):
            atoms_frac[i] = self.cart_to_frac(LAMMPS_structure.boxes[0], LAMMPS_structure.atoms[i])
            atoms_cart[i] = np.copy(LAMMPS_structure.atoms[i])
        #atoms_frac = np.copy(LAMMPS_structure.atoms)
        for i in range(len(atoms_frac)):
            ct = 0
            # We sort the typelist here to order the atoms in the output as desired
            for j in sorted_typelist:
                for at in range(len(atoms_frac[i])):
                    if LAMMPS_structure.id_lammps[i, at] == j:
                        sorted_atoms_frac[i, ct] = atoms_frac[i, at]
                        sorted_atoms_cart[i, ct] = atoms_cart[i, at]
                        ct += 1
        # We print a poscar here for free :)
        VASP_struct = VASP()
        VASP_struct.POSCAR.generate_poscar(B=np.array(LAMMPS_structure.boxes[0]), alat=1.0, names=sorted_namelist, mult=sorted_num_atoms, atoms_frac = sorted_atoms_frac[0])
        VASP_struct.POSCAR.write_frac(filename='POSCAR')

        self.print_sequential_MD_file(output_xdat, datalist=sorted_atoms_frac, headerMD=header + '\n')

    @calculate_time
    def vasp_to_xyz(self, VASP_structure):
        xyz_structure = xyz()
        # For now we only copy POSCARS sad :(
        structure = VASP_structure.POSCAR
        xyz_structure.title = structure.title + ''
        xyz_structure.B = np.copy(structure.B)
        xyz_structure.a_lat = structure.a_lat + 0.0
        xyz_structure.names = structure.names
        xyz_structure.multiplicity = np.copy(structure.multiplicity)
        xyz_structure.namelist = structure.namelist
        xyz_structure.atom_id = np.copy(structure.atom_id)
        xyz_structure.at_frac = np.copy(structure.at_frac)
        xyz_structure.at_cart = np.copy(structure.at_cart)

        self.type='xyz'
        return xyz_structure

    def vasp_to_scaleup(self, structure):
        pass
    
    def vasp_to_cif(self, structure):
        pass

    def vasp_to_res(self, structure):
        pass

    def scaleup_to_vasp(self, SCALEUP_structure):
        new_VASP = VASP()
        new_POSCAR = new_VASP.POSCAR

        new_POSCAR.title = SCALEUP_structure.title + ''
        new_POSCAR.B = np.copy(SCALEUP_structure.B)
        new_POSCAR.a_lat = SCALEUP_structure.a_lat + 0.0
        new_POSCAR.names = SCALEUP_structure.names
        new_POSCAR.lattype = 'Direct\n'
        new_POSCAR.multiplicity = np.copy(SCALEUP_structure.multiplicity)
        new_POSCAR.namelist = []

        # Atoms need to be reordered.
        ct = 0
        new_POSCAR.at_frac = np.zeros([len(SCALEUP_structure.at_frac), 3], dtype=float)
        new_POSCAR.at_cart = np.zeros([len(SCALEUP_structure.at_cart), 3], dtype=float)
        new_POSCAR.atom_id = np.zeros([len(SCALEUP_structure.atom_id)], dtype=float)
        for i in range(SCALEUP_structure.nspecies):
            for at in range(SCALEUP_structure.nat):
                if(SCALEUP_structure.at_type[at] == i):
                    new_POSCAR.at_frac[ct] = np.copy(SCALEUP_structure.at_frac[at])
                    new_POSCAR.at_cart[ct] = np.copy(SCALEUP_structure.at_cart[at])
                    new_POSCAR.atom_id[ct] = np.copy(SCALEUP_structure.atom_id[at])
                    new_POSCAR.namelist.append(SCALEUP_structure.namelist[at])
                    ct += 1
        new_POSCAR.type = "VASP"

        return new_VASP

    def GROMACStraj_to_vasp(self, GROMACS_structure, epoch=0):
        # Timestep here indicates the instance of GROMACS trajectory
        # that needs to be extracted
        new_VASP = VASP()
        new_POSCAR = new_VASP.POSCAR

        new_POSCAR.title = 'Gen. by Toolkit\n'
        new_POSCAR.B = np.copy(GROMACS_structure.boxes[epoch])
        new_POSCAR.a_lat = 1.0

        at_species = []
        spec_num = {}
        for spec in GROMACS_structure.atom_names[epoch]:
            if not spec in at_species:
                at_species.append(spec)
                spec_num[spec] = 1
            else:
                spec_num[spec] += 1
        new_POSCAR.names = np.array(at_species, dtype=str)
        new_POSCAR.lattype = 'Direct\n'
        new_POSCAR.multiplicity = np.array([spec_num[spec] for spec in at_species], dtype=int)
        new_POSCAR.at_frac = np.zeros([len(GROMACS_structure.atposfrac[epoch]), 3], dtype=float)
        new_POSCAR.at_cart = np.zeros([len(GROMACS_structure.atposfrac[epoch]), 3], dtype=float)
        new_POSCAR.atom_id = np.zeros([len(GROMACS_structure.atposfrac[epoch])], dtype=float)
        new_POSCAR.namelist = []
        ct = 0
        for spec in at_species:
            for at in range(len(GROMACS_structure.atposfrac[-1])):
                if GROMACS_structure.atom_names[-1, at] == spec:
                    new_POSCAR.at_frac[ct] = np.copy(GROMACS_structure.atposfrac[epoch, at])
                    new_POSCAR.at_cart[ct] = np.copy(GROMACS_structure.atposcart[epoch, at])
                    new_POSCAR.atom_id[ct] = np.copy(GROMACS_structure.index[epoch, at])
                    new_POSCAR.namelist.append(spec)
                    ct += 1

        new_POSCAR.type = "VASP"

        return new_VASP


    def GROMACS_to_vasp(self, GROMACS_structure):
        new_VASP = VASP()
        new_POSCAR = new_VASP.POSCAR

        new_POSCAR.title = 'Gen. by Toolkit\n'
        new_POSCAR.B = np.copy(GROMACS_structure.boxes[-1])
        new_POSCAR.a_lat = 1.0
        
        at_species = []
        spec_num = {}
        for spec in GROMACS_structure.atom_names[-1]:
            if not spec in at_species:
                at_species.append(spec)
                spec_num[spec] = 1
            else:
                spec_num[spec] += 1
        new_POSCAR.names = np.array(at_species, dtype=str)
        new_POSCAR.lattype = 'Direct\n'
        new_POSCAR.multiplicity = np.array([spec_num[spec] for spec in at_species], dtype=int)
        new_POSCAR.at_frac = np.zeros([len(GROMACS_structure.atposfrac[-1]), 3], dtype=float)
        new_POSCAR.at_cart = np.zeros([len(GROMACS_structure.atposfrac[-1]), 3], dtype=float)
        new_POSCAR.atom_id = np.zeros([len(GROMACS_structure.atposfrac[-1])], dtype=float)
        new_POSCAR.namelist = []
        ct = 0
        for spec in at_species:
            for at in range(len(GROMACS_structure.atposfrac[-1])):
                if GROMACS_structure.atom_names[-1, at] == spec:
                    new_POSCAR.at_frac[ct] = np.copy(GROMACS_structure.atposfrac[-1, at])
                    new_POSCAR.at_cart[ct] = np.copy(GROMACS_structure.atposcart[-1, at])
                    new_POSCAR.atom_id[ct] = np.copy(GROMACS_structure.index[-1, at])
                    new_POSCAR.namelist.append(spec)
                    ct += 1

        new_POSCAR.type = "VASP"

        return new_VASP

    def res_to_vasp(self, structure):
        pass

    def cif_to_vasp(self, structure):
        pass
    @calculate_time
    def transform_type(self, structure, target_type=None):
        # Controll the trolls
        if not target_type: print('usage: tools.transform_type(structure, "target_type")'); return -1;
        original_type = structure.type
        if original_type == target_type: print('Target type and original type are the same! Yo!'); return -1;

        # Target decider
        # Rule is: all goes to VASP if we can,
        # and then from VASP we go anywhere, limit is the sky.
        # Hence, first, type->VASP:
        if original_type != 'VASP':
            if original_type == 'scaleup':
                vasp_structure = self.scaleup_to_vasp(structure)
            if original_type == 'res':
                vasp_structure = self.res_to_vasp(structure)
            if original_type == 'cif':
                vasp_structure = self.cif_to_vasp(structure)
            if original_type == 'GROMACS':
                vasp_structure = self.GROMACS_to_vasp(structure)
        else:
            vasp_structure = structure
        
        # Thence, second, VASP->target_type:
        if target_type == 'xyz':
            new_structure = self.vasp_to_xyz(vasp_structure)
        if target_type == 'VASP':
            new_structure = vasp_structure
        if target_type == 'scaleup':
            new_structure = self.vasp_to_scaleup(vasp_structure)
        if target_type == 'cif': 
            new_structure = self.vasp_to_cif(vasp_structure)
        if target_type == 'res':
            new_structure = self.vasp_to_res(vasp_structure)

        return new_structure

