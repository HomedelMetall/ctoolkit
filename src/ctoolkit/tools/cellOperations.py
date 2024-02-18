from ctoolkit.global_vars.ext_libs import *

class cellOperations:
    def __init__(self):
        pass

    @calculate_time
    def func3(self):
        time.sleep(3)   

    @calculate_time
    def frac_to_cart(self, B, atom_set):
        new_atom_set = np.dot(B, atom_set.T).T
        
        return new_atom_set

    @calculate_time
    def cart_to_frac(self, B, atom_set):
        Binv = np.linalg.inv(B)
        new_atom_set = np.dot(Binv, atom_set.T).T

        return new_atom_set

    def rigid_displacement(self, vector, structure):
        new_cart = np.zeros([len(structure.at_cart), 3], dtype=float)
        new_frac = np.zeros([len(structure.at_frac), 3], dtype=float)
        vector = np.array(vector, dtype=float)

        disp = np.zeros([3], dtype=float)
        for i in range(3):
            disp += vector[i]*structure.B[i,:]

        for atom in range(len(structure.at_cart)):
            new_cart[atom] = np.copy(structure.at_cart[atom] + disp)
            
        structure.at_cart = np.copy(new_cart)
        structure.recalculate_poscar()

    def apply_rotation(self, vector, angle, axis, origin=None):
        # The idea: 
        # We want to rotate a group of atoms. The positions
        # of the atoms are given as vectors from the crystallographic
        # origin of coordinates. Thus what we want to do is:
        #   move origin -> apply rotation along given axis -> go back to crystal. origin

        # We use the tensorial form of the rotation operator
        axis = np.array(axis)
        axis = axis/(np.sqrt(np.dot(axis, axis)))
        cross = np.zeros([3,3], dtype=float)
        cart = np.array([[1.0, .0, .0], [.0, 1.0, .0], [.0, .0, 1.0]], dtype=float)
        for i in range(3):
            c1 = np.cross(axis, cart[i])
            cross += np.outer(c1, cart[i])
        R = np.cos(angle)*np.eye(3,3) + (np.sin(angle))*(cross) + (1-np.cos(angle))*(np.outer(axis, axis))

        # Now we move the origin, compute, and go back to the crystal. origin. All packed.
        if origin is None:
            origin = np.zeros([3], dtype=float)
        
        return np.dot(R, vector - origin) + origin
        
    def rigid_rotation(self, vector_list, angle, axis, origin=None):
        if len(vector_list.shape) < 2:
            return self.apply_rotation(vector_list, angle, axis, origin)
        else:
            veclist = []
            for i in range(len(vector_list)):
                veclist.append(self.apply_rotation(vector_list[i], angle, axis, origin))
            return np.array(veclist, dtype=float)

    def new_supercell_VASP(self, N, VASP_struct):
        new_VASP = self.copy_structure(VASP_struct)
        # Objects can be different...

        new_structure = new_VASP.POSCAR
        new_structure.at_frac = np.zeros([np.prod(N)*len(VASP_struct.POSCAR.at_frac), 3])
        ct_atom = 0
        ct2 = 0
        new_structure.namelist = ['']
        new_structure.atom_id =  []
        for specie_id, specie in enumerate(VASP_struct.POSCAR.names):
            new_structure.multiplicity[specie_id] = VASP_struct.POSCAR.multiplicity[specie_id]*np.prod(N)
            for atom in range(VASP_struct.POSCAR.multiplicity[specie_id]):
            #for nx, ny, nz in [(nx, ny, nz) for nx in range(N[0]) for ny in range(N[1]) for nz in range(N[2])]:
                #for atom in range(structure.multiplicity[specie_id]):
                for nx, ny, nz in [(nx, ny, nz) for nx in range(N[0]) for ny in range(N[1]) for nz in range(N[2])]: #This inverts how the atoms are printed
                    Id = np.eye(3)
                    cell_vector = nx*Id[0] + ny*Id[1] + nz*Id[2]

                    # Frac coordinates need to be reescaled to the new supercell
                    temp_vec = np.copy(VASP_struct.POSCAR.at_frac[ct2] + cell_vector)
                    for i in range(3):
                        new_structure.at_frac[ct_atom, i] = temp_vec[i]/N[i]
                    new_structure.namelist.append(VASP_struct.POSCAR.namelist[ct2])
                    new_structure.atom_id.append(VASP_struct.POSCAR.atom_id[ct2])
                    ct_atom += 1
                ct2 += 1

        for i in range(3):
            new_structure.B[i] = np.copy(VASP_struct.POSCAR.B[i]*N[i])

        new_structure.at_cart = np.copy(self.frac_to_cart(new_structure.B, new_structure.at_frac))

        mult = float(len(VASP_struct.POSCAR.at_frac))/float(len(new_structure.at_frac))
        new_structure.volume *= mult
        #new_structure.energy *= mult

        return new_structure

    def new_supercell(self, N, structure):
        # We select and divide into smaller functions as a sort of override
        if structure.type == 'VASP':
            new_structure = self.new_supercell_VASP(N, structure)

        return new_structure

    # Add atoms given a vector
    def spawn_cell_frac(self, vector, structure):
        new_frac = np.copy(structure.at_frac)
        for i in range(len(new_frac)):
            new_frac[i] = structure.at_frac[i] + np.array(vector, dtype=float)

        return np.copy(new_frac)

    def spawn_cell_invfrac(self, vector, structure):
        new_frac = np.copy(structure.at_frac)
        tol = 0.02
        for i in range(len(new_frac)):
            """
            temp_vec = -structure.at_frac[i] + np.array([1,1,1], dtype=float)
            # Correct frac on the fly
            for i in range(3):
                if(np.abs(temp_vec[i]-1.0) < tol):
                    temp_vec[i] -= 0.0#1.0

            #new_frac[i] = temp_vec + np.array(vector, dtype=float)
            """

            # We go to [-1,-1,-1] and then invert!
            temp_vec = -(structure.at_frac[i] + np.array([-1,-1,-1], dtype=float))
            # Correct frac on the fly
            for j in range(3):
                if(np.abs(temp_vec[j]-1.0) < tol):
                    temp_vec[j] -= 1.0

            new_frac[i] = temp_vec + np.array(vector, dtype=float)
            #new_frac[i] = np.array(vector, dtype=float)-structure.at_frac[i] + np.array([-1,-1,-1], dtype=float)

        return np.copy(new_frac)

    def centrosymmetric_supercell_VASP(self, VASP_structure):
        new_VASP = self.copy_structure(VASP_structure)
        new_structure = new_VASP.POSCAR

        # Centrosymmetry requires an atom well located for the origin.
        # Thus we need to shift to set the zero (for now experimental)
        shifted_structure = self.copy_structure(structure).POSCAR
        for i in range(len(shifted_structure.at_frac)):
            shifted_structure.at_frac[i] -= structure.at_frac[0]

        cells = []
        """
        for i, j in [(i,j) for i in range(2) for j in range(2)]:
            cells.append(self.spawn_cell_frac([i,j,0], shifted_structure))
        for i, j in [(i,j) for i in range(2) for j in range(2)]:
            cells.append(self.spawn_cell_invfrac([i,j,1], shifted_structure))
        """
        # Test: only 2 cells
        cells.append(self.spawn_cell_frac([0,0,0], shifted_structure))
        cells.append(self.spawn_cell_invfrac([0,0,1], shifted_structure))
        N = [1,1,2]
        ct = 0
        ct2 = 0
        new_structure.at_frac = np.zeros([np.prod(N)*len(structure.at_frac), 3])
        new_structure.namelist = ['']
        new_structure.atom_id = []
        for specie_id, specie in enumerate(structure.names):
            new_structure.multiplicity[specie_id] = structure.multiplicity[specie_id]*np.prod(N)
            for atom in range(structure.multiplicity[specie_id]):
                for cell in cells:
                    # Frac coordinates need to be reescaled to the new supercell
                    temp_vec = np.copy(cell[ct2])
                    for i in range(3):
                        new_structure.at_frac[ct, i] = temp_vec[i]/N[i]
                    new_structure.atom_id.append(structure.atom_id[ct2])
                    new_structure.namelist.append(structure.namelist[ct2])
                    ct += 1
                ct2 += 1

        for i in range(3):
            new_structure.B[i] = np.copy(structure.B[i]*N[i])

        new_structure.at_cart = np.copy(self.frac_to_cart(new_structure.B, new_structure.at_frac))

        mult = float(len(VASP_struct.POSCAR.at_frac))/float(len(new_structure.at_frac))
        new_structure.volume *= mult
        new_structure.energy *= mult

        return new_structure

    def centrosymmetric_supercell(self, structure):
        # We select and divide into smaller functions as a sort of override
        if structure.type == 'VASP':
            new_structure = self.new_supercell_VASP(N, structure)

        return new_structure

    def interpolate_structures_VASP(self, struct1, struct2, num_images):
        structures = []
        structures.append(struct1)
        for i in range(1, num_images+1):
            # Load structure + cleanup of frac coordinates
            new_VASP = self.copy_structure(struct1)
            new_structure = new_VASP.POSCAR
            new_structure.frac_image_sniper()
            # Linear interpolation of cell & frac & cart

            new_structure.B = struct1.POSCAR.B + float(i)/(num_images+2)*(struct2.POSCAR.B-struct1.POSCAR.B)
            new_structure.at_frac = struct1.POSCAR.at_frac + float(i)/(num_images+2)*(struct2.POSCAR.at_frac-struct1.POSCAR.at_frac)
            # Could put a control for running images here
            new_structure.at_cart = self.frac_to_cart(new_structure.B, new_structure.at_frac)#struct1.at_cart + float(i/(num_images+2))*(struct2.at_cart-struct1.at_cart)
            # Probably useless or troublesome:
            new_structure.frac_image_sniper()
            structures.append(new_VASP)

        structures.append(struct2)

        return structures

    def interpolate_structures(self, struct1, struct2, num_images):
        # We select and divide into smaller functions as a sort of override
        if struct1.type != struct2.type:
            print('Error! Structures in different format. Transform them first!')
            sys.exit()
            # Attempt conversion
            #print('Attempting to unify types...')
            #struct1 = self.transform_type(struct1, struct2.type)
            #print("Successfully transformed type of structure 1 (%s) into %s!" % (struct1.type, struct2.type))

        if struct1.type == 'VASP':
           structures = self.interpolate_structures_VASP(struct1, struct2, num_images)

        return structures

    @calculate_time
    # Tolerance is given in percentage
    def find_firstneighbors(self, atom, satellites, basis, tol=20):
        d = np.zeros([len(satellites)], dtype=float)
        #tol = 20 # in percentage
        for i in range(len(satellites)):
            best_transform, transform_key = self.findBestPairPeriodic(atom, satellites[i], basis)
            dist = np.sqrt(np.dot(atom-best_transform, atom-best_transform))
            # We exclude the "self" atom
            if(dist > 0.0):
                d[i] = dist
            else:
                d[i] = 1E5
        min_dist = np.min(d)
        list_neighbors = []
        for i in range(len(satellites)):
            if d[i] <= min_dist*(1+tol/100):
                #print(i, d[i], min_dist*(1+tol/100))
                list_neighbors.append(i)

        #print(list_neighbors, min_dist, np.sort(d))
        return list_neighbors

    def find_full_molecules(self, struct, center, satellite, num_neighbors=0):
        center_atoms = struct.filter_atoms(center)
        sat_atoms = struct.filter_atoms(satellite)
        all_atoms = np.copy(struct.at_cart)

        center_atoms_list = struct.filter_atoms_list(center)
        sat_atoms_list = struct.filter_atoms_list(satellite)

        groups = []
        nice_molecules = []
        for i in range(len(center_atoms)):
            nice_molecule = []
            sat_atom_neighbors = []
            list_neighbors = self.find_firstneighbors(center_atoms[i], sat_atoms, struct.B)

            for j in list_neighbors:
                sat_atom_neighbors.append(sat_atoms_list[j])

            if num_neighbors==1:
                # center atom new first neighbors
                cfn_bulk = self.find_firstneighbors(center_atoms[i], all_atoms, struct.B)
                cfn_list = []
                for cfn in cfn_bulk:
                    if ((cfn not in center_atoms_list) and
                        (cfn not in sat_atom_neighbors)):
                        cfn_list.append(cfn)

                # find the first neighbors of every satellite
                sfn_list_col = []
                for j in sat_atom_neighbors:
                    sfn_bulk = self.find_firstneighbors(struct.at_cart[j], all_atoms, struct.B)
                    sfn_list = []
                    for sfn in sfn_bulk:
                        if ((sfn not in center_atoms_list) and 
                            (sfn not in sat_atom_neighbors) and
                            (sfn not in cfn_list)): 
                            sfn_list.append(sfn)

                    sfn_list_col.append(sfn_list)
                 
                groups.append([center_atoms_list[i], sat_atom_neighbors, cfn_list, sfn_list_col])
                nice_molecule.append(center_atoms_list[i])
                for sat in sat_atom_neighbors: nice_molecule.append(sat)
                for cfn in cfn_list: nice_molecule.append(cfn)
                for sfn in sfn_list: nice_molecule.append(sfn)
                nice_molecules.append(nice_molecule)
            else:
                groups.append([center_atoms_list[i], sat_atom_neighbors])

            sys.stdout.write("%3.3f%%\r" % (100*(i+1)/len(center_atoms)))
            sys.stdout.flush()

        return groups, nice_molecules

    # Find groups of atoms between the CENTER
    # of the interaction and the SATELLITE possible
    # atoms belonging to a particular species
    @calculate_time
    def find_groups(self, struct, center, satellite, excludeSatList=[], tol=20):
        center_atoms = struct.filter_atoms(center)
        sat_atoms = struct.filter_atoms(satellite, excludeSatList=excludeSatList)
        center_atoms_list = struct.filter_atoms_list(center)
        sat_atoms_list = struct.filter_atoms_list(satellite, excludeSatList=excludeSatList)
        #print(len(sat_atoms), len(sat_atoms_list))
        groups = []
        for i in range(len(center_atoms)):
            sat_atom_neighbors = []
            list_neighbors = self.find_firstneighbors(center_atoms[i], sat_atoms, struct.B, tol=tol)
            for j in list_neighbors:
                sat_atom_neighbors.append(sat_atoms_list[j])
            groups.append([center_atoms_list[i], sat_atom_neighbors])
            
            #sys.stdout.write("%3.3f%%\r" % (100*(i+1)/len(center_atoms)))
            #sys.stdout.flush()

        return groups

    # Sometimes atoms go out of the simulation box. Bound box. Whatever.
    # And sometimes we want them out to understand the structure.
    # So I add here a wrapper for:
    #   1) Application of periodic boundary conditions.
    #   2) Unbound atoms.
    # The code here works in fractional coordinates. 
    # You can provide a structure, or a set of atoms,
    # or a cartesian set of atoms + basis...
    # Or at least, it will do that. For now, at_frac.
    @calculate_time
    def handle_periodic_boundaries(self, at_frac, mode='bound'):
        if mode != 'bound' or mode != 'unbound':
            print("PBC handler: using mode 'bound' by default (input not recognized)")
            mode = 'bound'
            #CONTINUE_HERE

    # This function searches for a periodic image of periodic_atom
    # that is the nearest to fixed_atom. Returns the best periodic
    # image found, and the key [k_x, k_y, k_z] necessary to transform it.
    @calculate_time
    def findBestPairPeriodic_test(self, fixed_atom, periodic_atom, basis):
        mindist = 1E3
        best_transform = np.copy(periodic_atom)
        best_transform_key = np.array([0,0,0], dtype=float)
        # We search 1 cell around the boundary
        # Maybe we can search in a decomposed way i>j>k?????
        for axis in range(3):
            transform_key = np.copy(best_transform_key)
            for i in range(-1, 2):
                transform_key[axis] = i
                transformed_atom = periodic_atom + basis[0]*transform_key[0] + basis[1]*transform_key[1] + basis[2]*transform_key[2]
                vec = fixed_atom-transformed_atom
                dist = np.sqrt(np.dot(vec, vec))
                if (dist-mindist) < 0.0:
                    mindist = dist + 0.0
                    best_transform = np.copy(transformed_atom)
                    best_transform_key[axis] = i

        return best_transform, best_transform_key

    @calculate_time
    def findBestPairPeriodic_GPT(self, fixed_atom, periodic_atom, basis):
        translations = np.array([i*basis[0] + j*basis[1] + k*basis[2] for i in range(-1,2) for j in range(-1,2) for k in range(-1,2)])
        #transformed_atoms = periodic_atom + translations
        #print(translations)
        transformed_atoms = np.array([periodic_atom + translations[i] for i in range(len(translations))])
        #print(transformed_atoms)
        #vec = fixed_atom - transformed_atoms
        vec = np.array([fixed_atom + transformed_atoms[i] for i in range(len(transformed_atoms))])
        dist = np.linalg.norm(vec, axis=1)
        min_index = np.argmin(dist)
        #print(np.min(dist), min_index)
        return transformed_atoms[min_index], [min_index//9-1, (min_index//3)%3-1, min_index%3-1]

    # This function searches for a periodic image of periodic_atom
    # that is the nearest to fixed_atom. Returns the best periodic
    # image found, and the key [k_x, k_y, k_z] necessary to transform it.
    @calculate_time
    def findBestPairPeriodic(self, fixed_atom, periodic_atom, basis):
        mindist = 1E3
        best_transform = np.copy(periodic_atom)
        transform_key = [0,0,0]
        # We search 1 cell around the boundar
        """
        bi, bj, bk = 0, 0, 0
        for i in range(-1, 2):
            transformed_atom = periodic_atom + basis[0]*float(i) + basis[1]*float(bj) + basis[2]*float(bk)
            vec = fixed_atom-transformed_atom
            dist = np.sqrt(np.dot(vec, vec))
            try:
                if (dist-mindist) < 0.0:
                    mindist = dist
                    best_transform = np.copy(transformed_atom)
                    bi = i
            except:
                print("IT BROKE! findBestPeriodic doesn't work well!", dist, mindist, vec)
                sys.exit()

        for j in range(-1, 2):
            transformed_atom = periodic_atom + basis[0]*float(bi) + basis[1]*float(j) + basis[2]*float(bk)
            vec = fixed_atom-transformed_atom
            dist = np.sqrt(np.dot(vec, vec))
            try:
                if (dist-mindist) < 0.0:
                    mindist = dist
                    best_transform = np.copy(transformed_atom)
                    bj = j
            except:
                print("IT BROKE! findBestPeriodic doesn't work well!", dist, mindist, vec)
                sys.exit()

        for k in range(-1, 2):
            transformed_atom = periodic_atom + basis[0]*float(bi) + basis[1]*float(bj) + basis[2]*float(k)
            vec = fixed_atom-transformed_atom
            dist = np.sqrt(np.dot(vec, vec))
            try:
                if (dist-mindist) < 0.0:
                    mindist = dist
                    best_transform = np.copy(transformed_atom)
                    bk = k
            except:
                print("IT BROKE! findBestPeriodic doesn't work well!", dist, mindist, vec)
                sys.exit()

        transform_key = [bi, bj, bk]
        """
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    transformed_atom = periodic_atom + basis[0]*float(i) + basis[1]*float(j) + basis[2]*float(k)                
                    vec = fixed_atom-transformed_atom
                    dist = np.sqrt(np.dot(vec, vec))
                    #print(transformed_atom, dist, mindist, best_transform)
                    try:
                        if (dist-mindist) < 0.0:
                            mindist = dist
                            best_transform = np.copy(transformed_atom)
                            transform_key = [i, j, k]
                    except:
                        print("IT BROKE! findBestPeriodic doesn't work well!", dist, mindist, vec)
                        sys.exit()
        
        return best_transform, transform_key

    @calculate_time
    def atom_transform(self, atom, basis, key):
        return atom + basis[0]*float(key[0]) + basis[1]*float(key[1]) + basis[2]*float(key[2])

    # Function to find pairs of atoms. Input can be in standard .at_frac
    # Returns the pair[atom1]=atom2 and the keyring [k_x, k_y, k_z] for the periodic images.
    @calculate_time
    def find_pairs(self, atoms1, atoms2, box2):
        pair = np.zeros([len(atoms1)], dtype=int)
        pair[:] -= 1
        mindist_pair = np.zeros([len(atoms1)], dtype=float)
        mindist_pair[:] -= 1
        keyring = []
        listat2 = list(range(len(atoms2)))
        for i in range(len(atoms1)):
            threshold = 1E3 # This is just an init parameter
            for j in listat2:
                atoms2periodic, transform_key = self.findBestPairPeriodic_test(atoms1[i], atoms2[j], box2)
                dvec = np.copy(atoms1[i]-atoms2periodic)
                dist = np.sqrt(np.dot(dvec, dvec))
                # This is unclean but efficient. Basically, checks that the new transformation
                # provides a smaller distance between pairs of atoms. Probably can be cleaned.
                if (dist-threshold) < 0.0:
                    threshold = dist + 0.0
                    pair[i] = j
                    best_transform_key = transform_key
            if pair[i] in listat2: listat2.remove(pair[i])
            keyring.append(best_transform_key)
            mindist_pair[i] = threshold + 0.0

        for i, a2 in enumerate(pair):
            # This warning should be set to an actual threshold parameter
            if mindist_pair[i] > 2.5 and a2 > -1: print("WARNING with pair %d %d min distance: %.5f" % (i, a2, mindist_pair[i]))
            #print("Pair %d %d: %.5f" % (i, a2, mindist_pair[i]), atoms1[i], atoms2[a2], keyring[i]) # DEBUG

        return pair, keyring

    # A quicker function to find keys of known pairs of atoms
    @calculate_time
    def refresh_key(self, atoms1, atoms2, box2, pairs):
        key = []
        for at1, at2 in enumerate(pairs):
            atoms2periodic, transform_key = self.findBestPairPeriodic_test(atoms1[at1], atoms2[at2], box2)
            key.append(transform_key)
        return key

    # Key is for the atom 2.
    @calculate_time
    def build_vectormap(self, atoms1, atoms2, box, pairs, initial_key):
        vectors = np.zeros([len(atoms1), 3], dtype=float)
        key = initial_key + .0
        d = []
        for i, j in enumerate(pairs):
            atoms2_pos_correction = key[i][0]*box[0] + key[i][1]*box[1] + key[i][2]*box[2]
            vec = atoms1[i] - (atoms2[j] + atoms2_pos_correction)
            vectors[i] = np.copy(vec)
            d.append(np.sqrt(np.dot(vec,vec)))
            #print(d[-1]) 
        # The tolerance should probably be a parameter
        if (np.max(np.array(d, dtype=float)) > 2.5):
            key = self.refresh_key(atoms1, atoms2, box, pairs)
            d = []
            for i, j in enumerate(pairs):
                atoms2_pos_correction = key[i][0]*box[0] + key[i][1]*box[1] + key[i][2]*box[2]
                vec = atoms1[i] - (atoms2[j] + atoms2_pos_correction)
                vectors[i] = np.copy(vec)
                dist = np.sqrt(np.dot(vec, vec))
                d.append(dist)
        
        # Just a check
        if (np.max(np.array(d, dtype=float)) > 2.5):
            print("Max distance between atoms is larger than 2.5")
            sys.exit()

        return vectors, key

    def compute_angles_MD_GPT(self, vectors, vec_of_vecs=False):
        num_steps, num_vecs = vectors.shape[:2]
        e_x = np.array([1, 0, 0], dtype=np.float32)
        e_z = np.array([0, 0, 1], dtype=np.float32)
        v_pbc = vectors - np.einsum('ij,j->ij', vectors, e_x) * e_x
        v_pab = vectors - np.einsum('ij,j->ij', vectors, e_z) * e_z
        norm_vpbc = np.linalg.norm(v_pbc, axis=-1)
        norm_vpab = np.linalg.norm(v_pab, axis=-1)
        phi = np.arctan2(np.cross(v_pab, e_z), np.dot(v_pab, e_x))
        theta = np.arccos(np.einsum('ij,j->i', vectors, e_z) / norm_vpbc)
        ctheta = np.cos(theta)
        if vec_of_vecs:
            angles = np.concatenate((phi.reshape(-1,1), theta.reshape(-1,1)), axis=1)
            return angles.reshape(num_steps, num_vecs, 2)
        else:
            return phi, theta, ctheta

    @calculate_time
    def compute_angles_MD(self, vectors, vec_of_vecs=False):
        num_steps = len(vectors)
        num_vecs = len(vectors[0])
        theta = np.zeros([num_steps*num_vecs], dtype=float)
        phi = np.zeros([num_steps*num_vecs], dtype=float)
        ctheta = np.zeros([num_steps*num_vecs], dtype=float)

        e_x = np.eye(3, dtype=float)[0,:]
        e_z = np.eye(3, dtype=float)[2,:]
        angles = []
        for i in range(num_steps):
            angles_inner = []
            for j in range(num_vecs):
                vec = np.copy(vectors[i,j])

                v_pbc = vec - np.dot(vec, e_x)*e_x
                v_pab = vec - np.dot(vec, e_z)*e_z
                norm_vpbc = np.sqrt(np.dot(v_pbc, v_pbc))
                norm_vpab = np.sqrt(np.dot(v_pab, v_pab))
                # Angles = [\phi, \theta]
                phi[i*num_vecs+j] = self.calculate_angle_full(e_x, v_pab/norm_vpab, e_z) + 0.0
                theta[i*num_vecs+j] = self.calculate_angle(e_z, vec)#angle(e_z, v_pbc/norm_vpbc, e_x) + 0.0
                ctheta[i*num_vecs+j] = np.cos(theta[i*num_vecs+j])
                angles_inner.append([phi[i*num_vecs+j], theta[i*num_vecs+j]])
            angles.append(angles_inner)

        angles = np.array(angles)

        if vec_of_vecs == False:
            return phi, theta, ctheta
        else:
            return angles

