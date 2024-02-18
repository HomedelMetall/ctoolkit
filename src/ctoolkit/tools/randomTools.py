from ctoolkit.global_vars.ext_libs import *
from ctoolkit.global_vars.tk_lib_LAMMPS import *
from ctoolkit.global_vars.tk_lib_VASP import *
class randomTools():
    def __init__(self):
        pass

    def cat_lammps_trajectories(self, filenames):
        output_xdat = 'XDATCAR'
        LAMMPS_structs = []
        total_nat = 0
        for name in filenames:
            Lstruct = LAMMPS()
            Lstruct.read_snapshots(name, style='typepos')
            total_nat += len(Lstruct.id_lammps[0])
            LAMMPS_structs.append(Lstruct)

        LAMMPS_struct_cat = LAMMPS()
        LAMMPS_struct_cat.id_lammps = np.zeros([len(LAMMPS_structs[0].id_lammps),
                                                    total_nat], dtype=int)
        LAMMPS_struct_cat.boxes = np.copy(LAMMPS_structs[0].boxes)
        LAMMPS_struct_cat.atoms = np.zeros([len(LAMMPS_structs[0].atoms),
                                                    total_nat, 3], dtype=float)

        init = 0
        for struct in LAMMPS_structs:
            LAMMPS_struct_cat.id_lammps[:, init:init+len(struct.id_lammps[0])] = np.copy(struct.id_lammps)
            LAMMPS_struct_cat.atoms[:, init:init+len(struct.atoms[0]), :] = np.copy(struct.atoms)
            init += len(struct.id_lammps[0])

        self.lammps_to_xdatcar(LAMMPS_struct_cat, output_xdat)

    def group_coherence_test(self, group):
        chk = len(group[0][1])
        for g in group:
            if len(g[1]) != chk:
                print("Not all centers have the same number of group satellites! Check tolerances...")
                return -1

        return 0

    def molecule_finder(self, species):
        VASPstruct = VASP()
        VASPstruct.POSCAR.read_poscar('POSCAR')
        VASP_atoms = np.copy(VASPstruct.POSCAR.at_cart)
        #MD_structures = tools.loadVASPMD('.', 'XDATCAR')
        print('Looking for groups...')
        groups = []
        for ids, specie1 in enumerate(species):
            for specie2 in species[ids:]: 
                group = self.find_groups(VASPstruct.POSCAR, specie1, specie2)
                print('%s-%s: %d groups of %d atoms' % (specie1, specie2, len(group), 1+len(group[0][1])))
                result = group_coherence_test(group)

                groups.append(group)

    @calculate_time
    def readXDATCAR(self, filename):
        with open(filename, 'r') as file:
            header_lines = 7
            B = np.zeros([3,3], dtype=float)
            step_header_lines = 1
            lines = file.readlines()
            for i in range(3):
              B[i,:] = np.array(lines[2+i].split()[:],dtype=float)
            numatoms = np.sum(np.array(lines[6].split(), dtype=int))
            numsnapshots = int((len(lines)-header_lines)/(step_header_lines + numatoms))
            print(numsnapshots) 
            at_frac = np.zeros([numsnapshots, numatoms, 3], dtype=float)
            for i in range(numsnapshots):
                pos_file = header_lines + step_header_lines + (step_header_lines+numatoms)*i
                at_frac[i] = np.array([[float(val) for val in line.split()] for line in lines[pos_file:pos_file+numatoms]])
        print(at_frac.shape)
        return B, at_frac

    @calculate_time
    def linearfit(self, x, y):
        # Solving for y = mx + c
        A = np.vstack([x, np.ones(len(x))]).T
        m, c = np.linalg.lstsq(A, y, rcond=None)[0]
        return m, c

    @calculate_time
    def sorting_order(self, array, tol=1E-8):
        sorted_list = np.sort(array)
        s_order = []
        for i in range(len(sorted_list)):
            for j in range(len(array)):
                if (np.abs(array[j]-sorted_list[i]) < tol) :
                    s_order.append(j)
                    break

        return s_order

    @calculate_time
    def resort(self, array, sorting_order):
        new_array = np.copy(array)
        for i, sord in enumerate(sorting_order):
            new_array[i] = array[sord]
    
        return new_array

    @calculate_time
    def sortMatrixbyCols(self, A, axis, tol=1E-8):
        Anew = np.copy(np.array(A, dtype=float))
        sorted_axis = np.copy(np.sort(Anew[:, axis]))
        Asorted = np.zeros([len(Anew), len(Anew[0])], dtype=float)
        for i in range(len(sorted_axis)):
            for j in range(len(Anew[:, axis])):
                if(np.abs(sorted_axis[i] - Anew[j, axis]) < tol):
                    for l in range(len(Anew[j, :])):
                        Asorted[i,l] = Anew[j,l] + .0

        return Asorted

    # Num_info = # atoms, for example
    @calculate_time
    def read_sequential_file(self, filename, num_info, hdr_lines = 0, offset = 0):
        # Basically: a sequential file has N header lines + M atom lines for S snapshots, for example.
        # Hence, the number of lines L has to be
        #           L = (N+M)*S

        with open(filename, 'r') as file:
            lines = file.readlines()
            num_steps = int((len(lines)-offset)/(hdr_lines + num_info))
            res = (len(lines)-offset)%(hdr_lines + num_info)
            if res != 0:
                print("Hey! Check num_info, hdr_lines. The data you provide seems non commensurate!\n Residual (L-offset)%%(N+M) = %.3f" % (res))
            dim_data = len(lines[offset+hdr_lines].split())
            data = np.zeros([num_steps, num_info, dim_data], dtype=float)
            hdr = []
            for i in range(num_steps):
                sp = int(offset+i*(hdr_lines+num_info))
                data[i] = np.array([[float(val) for val in line.split()] for line in lines[sp+hdr_lines:sp+hdr_lines+num_info]])
                #for j in range(hdr_lines):
                #    hdr.append(lines[sp+j].split())
                #for j in range(num_info):
                #    data[i,j] = np.array(lines[sp+hdr_lines+j].split(), dtype=float)
                sys.stdout.write("%3.2f...\r" % (100*(i+1)/num_steps))
                sys.stdout.flush()

        # Handle hdr with care. It's a list of np.arrays, not an np.array itself!
        return np.copy(data)#hdr, np.copy(data)

    @calculate_time
    def read_2D_file(self, filename, offset=None):
        with open(filename, 'r') as file:
            lines = file.readlines()
            
            if(offset == None):
                offset = 0
                while('#' in lines[offset]):
                    offset += 1

            num_rows = len(lines)-offset
            num_cols = len(lines[offset].split())
            data = np.zeros([num_rows, num_cols], dtype=float)
            for i, line in enumerate(lines[offset:]):
                data[i] = np.array(line.split(), dtype=float)

        return np.copy(data)

    @calculate_time
    def print_sequential_MD_file(self, filename, datalist, headerMD='#\n'):
        numsteps = len(datalist)
        numcols = len(datalist[0,0])
        s = ''
        s += headerMD
        if not '\n' in headerMD:
            s += '\n'

        for step in range(numsteps):
            header_step = 'Direct configuration=     %d\n' % (step)
            datanow = []
            for i in range(numcols):
                datanow.append(datalist[step,:,i])

            snow = self.print_to_file(filename, datanow, header=None, tostring=True)

            s += header_step + snow

        with open(filename, 'w') as file:
            file.write(s)

    @calculate_time
    def print_to_file_GPT(self, filename, datalist, header='#First line is comment line\n', tostring=False, precision=8):
        # Check that the length of the file to print in rows is consistent for all data...
        length_check = np.array([data.shape[0] for data in datalist], dtype=int)
        if len(np.unique(length_check)) != 1:
            raise ValueError("Data length of rows inconsistent!")
        # Use numpy's savetxt() method to write the data to a file
        np.savetxt(filename, np.column_stack(datalist), header=header, fmt='%.*e\t' % precision)

    # A wrapper for file print functions. 
    # datalist = [np.array1, ...]
    @calculate_time
    def print_to_file(self, filename, datalist, header='#First line is comment line\n', tostring=False, precision=8):
        fopen = open(filename, 'w')
        numcols = len(datalist)
        multiplicity = []
        length_check = []
        for data in datalist:
            # Values
            if len(data.shape) == 0:
                multiplicity.append(1)
                length_check.append(len(data))

            # Columns
            elif len(data.shape) == 1:
                multiplicity.append(1)
                length_check.append(len(data))

            # Matrices
            elif len(data.shape) > 1:
                (dimrows, dimcols) = data.shape
                length_check.append(dimrows)
                multiplicity.append(dimcols)
            
        # Check that the length of the file to print in rows is consistent for all data...
        length_check = np.array(length_check, dtype=int)
        for i in range(len(length_check)):
            if length_check[0] != length_check[i]:
                print("Data length of rows inconsistent!")
                sys.exit()

        # Usual data->string transformation :)
        s = ''
        for row in range(length_check[0]):
            for i, data in enumerate(datalist):
                if multiplicity[i] == 1:
                    s += '%.*e\t' % (precision, data[row])
                else:
                    for mult in range(multiplicity[i]):
                        s += '%.*e\t' % (precision, data[row, mult])
            s += '\n'
        
        #print("I'm here")
        if(tostring == True):
            return s
        
        with open(filename, 'w') as file:
            if(header == None):
                file.write(s)
            else:
                if not '\n' in header:
                    header += '\n'
                file.write(header+s)

    @calculate_time
    def flatten_3dvecs(self, array):
        darray = np.zeros([3*len(array)], dtype=float)
        ct = 0
        for i in range(len(array)):
            for j in range(len(array[0])):
                darray[ct] = array[i, j]

        return np.copy(darray)

    @calculate_time
    def blow_3dvecs(self, array):
        darray = np.zeros([len(array)/3, 3], dtype=float)
        ct = 0
        for i in range(len(array)):
            for j in range(3):
                darray[i,j] = array[3*i+j]

        return np.copy(darray)       

    @calculate_time
    def vector_of_norms(self, array):
        darray = np.zeros([len(array), 3], dtype=float)
        for i in range(len(array)):
            darray[i, 0] = 1. / np.sqrt(np.dot(array[i], array[i]))
            darray[i, 1] = 1. / np.sqrt(np.dot(array[i], array[i]))
            darray[i, 2] = 1. / np.sqrt(np.dot(array[i], array[i]))

        return np.copy(darray)

    @calculate_time
    def floatrange(self, start, end, numdiv=100):
        vec = np.zeros([numdiv], dtype=float)
        for i in range(numdiv):
            vec[i] = start + float(i*(end-start)/numdiv)

        return vec

    @calculate_time
    def normalize_vector_set(self, vectorset):
        num_vecs = len(vectorset)
        dims = len(vectorset[0])
        norms = 1./np.array(list(map(np.sqrt, np.einsum('...j,...j', vectorset, vectorset))), dtype=float)
        newvec = np.copy(np.einsum('...i,...', vectorset, norms))
    
        return np.copy(newvec)

    @calculate_time
    def normalize_vectors_MD(self, vectors):
        num_steps = len(vectors)
        num_vecs = len(vectors[0])
        dims = len(vectors[0,0])

        newvecs = np.zeros([num_steps, num_vecs, dims], dtype=float)

        for i in range(num_steps):
            newvecs[i] = np.copy(self.normalize_vector_set(vectors[i]))

        return newvecs

    # Average a set of lists, possibly applying an offset to data :)
    @calculate_time
    def multiaverage(self, collection, offset):
        output = []
        for item in collection:
            output.append(np.average(item[offset:len(item)]))

        return output

    def print_timers(self):
        timers_dict['Total run'] = time.time() - timers_dict['Total run']
        print("Run timers summary. Total run: %d secs" % (timers_dict['Total run']))
        print("===================================")
        for timer_name in dict(sorted(timers_dict.items(), key=lambda x:x[1], reverse=True)):
            print("Function %s: %d secs" % (timer_name, timers_dict[timer_name]))
    
    def plot(self, a, xrange=None, yrange=None, filename=None):
        import matplotlib.pyplot as plt
        # Assumption: each element of a is a two-item list.
        for item in a:
            plt.plot(item[0], item[1])

        ax = plt.gca()
        if(xrange != None):
            ax.set_xlim([xrange[0],xrange[1]])
        if(yrange != None):
            ax.set_ylim([yrange[0], yrange[1]])
        if(filename != None):
            plt.savefig(filename, format='png')
        if(filename == None):
            plt.show()

    def movie_from_images(image_path_list, output='output.avi', video_length_secs=None):
        import cv2
        frames = []
        for p in image_path_list:
            frames.append(cv2.imread(p))

        if video_length_secs is None:
            video_length_secs = len(frames)

        fps = len(frames) / video_length_secs
        height, width, layers = frames[0].shape
        video = cv2.VideoWriter(output, codec=cv2.CVideoWriter_fourcc(*'DIVX'), frames_per_second = fps, frame_size = (width, height))
        [video.write(frame) for frame in frames]
        cv2.destroyAllWindows()
        video.release()
