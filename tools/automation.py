from toolkit.global_vars.ext_libs import *
from toolkit.global_vars.tk_lib_VASP import *
import shutil

# Automated utilities
class automation:
    def __init__(self):
        pass
    
    @calculate_time
    def func2(self):
        time.sleep(2)

    # Given a structure of Directories/Subdirectories, where
    # each has a tag (for example, N_kbar/MK_run), read
    # the contents of files named "filename" and returns:
    #   - names: value of the subdirectory name, assuming it's linked to a physical property and its an integer
    #   - dataset: sets of the data stored in the "filename" file
    #   - rootdirs: supertag directory names
    @calculate_time
    def read_file_study(self, filename, supertag='kbar', subtag='run'):
        dirlist = sorted(os.listdir("./"))
        dataset = []
        names = []
        rootdirs = []
        # If we create a decider here for no subtag case, we can turn the function
        # into a 1-dirlevel study processing tool too.
        for dirs in dirlist:
            if (supertag in dirs) and os.path.isdir(dirs) is True:
                subdirlist = sorted(os.listdir("./"+dirs+"/"))
                subdataset = []
                subnames = []
                for subdirs in subdirlist:
                    if (subtag in subdirs) and os.path.isdir(dirs+"/"+subdirs):
                        if os.path.isfile(dirs+'/'+subdirs+'/'+filename) is False:
                            continue
                        data = self.read_2D_file(dirs+'/'+subdirs+'/'+filename)
                        subdataset.append(data)
                        # Assuming the subtag is about integer temperatures...
                        T = int(re.findall(r'[\d.]+', subdirs)[0])
                        subnames.append(T)

                # Datalength consistency check
                flag = False
                # Just a patch that assumes that first calculation is fine
                incompatible_sets = []
                for i in range(len(subdataset)):
                    if(len(subdataset[i]) != len(subdataset[0])):
                        flag = True
                        incompatible_sets.append(i)
                if flag == True:
                    for i in incompatible_sets:
                        del subnames[i]
                        del subdataset[i]

                sord = self.sorting_order(subnames)
                ssubnames = self.resort(subnames, sord)
                ssubdataset = self.resort(subdataset, sord)

                names.append(ssubnames)
                dataset.append(np.array(ssubdataset))
                rootdirs.append(dirs)

        return names, dataset, rootdirs

    # A way to load XDATCAR info for constant lattice parameter simulations.
    @calculate_time
    def loadVASPMD(self, sim_directory='', filename = None):
        structures = []
        print('Automatic load of VASP Molecular Dynamics nvt...')
        
        # We asume (for simplicity of the case of implementation)
        # that the files are called "XDATCARN", where N is the number that
        # orders the files. So we figure this out first:
        num_items = len(os.listdir(sim_directory))
        if filename != None:
            filenames = [sim_directory + '/' + filename]
        else:
            filenames = []
            for i in range(num_items+1):
                if(self.fileExists(sim_directory + '/XDATCAR%d' % i) == 0):
                    filenames.append(sim_directory + '/XDATCAR%d' % i)

        # In order to load the structures, I believe that what is probably more memory costly
        # but simpler is to create one VASP structure per snapshot. It will be clearer, otherwise
        # we have to modify the main class VASP so that it can contain several snapshots.
        # Therefore, we read first the XDATCARs

        MD_traj_frac = []
        for filename in filenames:
            print("Reading %s..." % (filename))
            basis, at_frac = self.readXDATCAR(filename)
            for i in range(len(at_frac)):
                MD_traj_frac.append(np.copy(at_frac[i]))
        MD_traj_frac = np.copy(np.array(MD_traj_frac, dtype=float))

        # We now load the POSCAR of the folder over a dummy VASP structure
        # And we span the all the new structures as copies of it
        # In case things have been generated with PHONOPY, we also look for POSCARini
        
        if os.path.isfile(sim_directory + '/POSCAR'):
            POSCAR = sim_directory + '/POSCAR'
        else:
            POSCAR = sim_directory + '/POSCARini'
        
        vaspDUM = VASP()
        if self.fileExists(POSCAR) == 0: vaspDUM.POSCAR.read_poscar(POSCAR)
        #vaspDUM.POSCAR.generate_poscar(B=basis, alat=1.0, names=sorted_namelist, mult=sorted_num_atoms, atoms_frac=sorted_atoms_frac[0])
        MD_structures = []
        for i in range(len(MD_traj_frac)):
            new_struct = self.copy_structure(vaspDUM)
            MD_structures.append(new_struct)

        print("Generating %d structures..." % (len(MD_traj_frac)))
        for i, struct in enumerate(MD_structures):
            struct.POSCAR.at_cart =  self.frac_to_cart(struct.POSCAR.B, MD_traj_frac[i])
            struct.POSCAR.recalculate_poscar()
            sys.stdout.write("%3.2f%%\r" % (100*(i+1)/len(MD_traj_frac)))
            sys.stdout.flush()
        
        # This is cool, but not enough. There can't be a good way of setting up image sniper tools with a system
        # that can rotate.

        # Instead, we standarize the periodic images to be always "near" each other. No big jumps. I.e., if
        # an atom crosses more than half a unit cell in a whim, it's probably needing correction.
        print("Standarizing %d structures..." % (len(MD_traj_frac)))
        std_MD = self.standarizeMDVASP(MD_structures)

        return std_MD

    @calculate_time
    def standarizeMDVASP(self, structs):
        tol = 0.85
        for i in range(len(structs)-1):
            frac_disp = structs[i].POSCAR.at_frac - structs[i+1].POSCAR.at_frac
            structs[i+1].POSCAR.at_frac[frac_disp < -tol] -= 1.0
            structs[i+1].POSCAR.at_frac[frac_disp > tol] += 1.0
            structs[i+1].POSCAR.at_cart = self.frac_to_cart(structs[i+1].POSCAR.B, structs[i+1].POSCAR.at_frac)

        return structs

    # Automatic load of VASP structures. It also
    # detects and processes the properties that
    # can be extracted from it
    @calculate_time
    def loadVASPStructs(self, dirlist=None):
        structures = []
        print('Automatic load of VASP calculations...')
        # The user may input their own list of directories to load...
        if dirlist is None: VASPdirs = self.getCalculationDirs()
        else: VASPdirs = dirlist
        print('Loading...')
        for dirs in VASPdirs:
            vaspNow = VASP()
            POSCAR = dirs + '/POSCAR'
            CONTCAR = dirs + '/CONTCAR'
            OUTCAR = dirs + '/OUTCAR'
            if self.fileExists(POSCAR) == 0: vaspNow.POSCAR.read_poscar(POSCAR)
            if self.fileExists(CONTCAR) == 0: vaspNow.POSCAR.read_poscar(CONTCAR)
            # Load all properties available from OUTCAR
            props = ['BC', 'E', 'V', 'AllE', 'AllV']
            propsDic = {'BC' : 'Born effective charges', 'E' : 'Energy', 'V' : 'Volume', 'AllE' : 'Energy trajectory', 'AllV' : 'Volume trajectory'}
            if self.fileExists(OUTCAR) == 0:
                props_okay = []
                for prop in props:
                    try:
                        vaspNow.outcar.extract_data(OUTCAR, prop)
                        props_okay.append(prop)
                    except:
                        continue

            structures.append(vaspNow)
            if(vaspNow.getStatus()['POSCAR']): POSCAR_status='OK'
            else: POSCAR_status = 'not read'
            if(vaspNow.getStatus()['OUTCAR']): 
                OUTCAR_status='OK'
                OUTCAR_tags = ''
                for i, prop in enumerate(props_okay):
                    OUTCAR_tags += prop
                    if i < len(props_okay)-1:
                        OUTCAR_tags += ', '
                print('\t%10s loaded! POSCAR %s. OUTCAR %s. Read from OUTCAR: %s.' % (dirs, POSCAR_status, OUTCAR_status, OUTCAR_tags))
            else:
                print('\t%10s loaded! POSCAR %s. OUTCAR not read.' % (dirs, POSCAR_status))
        return structures

    # Only 1 variable studies for now
    @calculate_time
    def createStudyVASP(self, var, structures, dirheader, style, launch_command='echo pwd'):
        print('Preparing directories for study...')
        dirnames = self.spanDirs(var, dirheader, style)
        if self.copyBASEdir(dirnames) == -1:
            print("Please, setup correctly the environment and try again.")
            for dirname in dirnames:
                shutil.rmtree(dirname)
            sys.exit()
        
        print('INPUT_BASE copied')
        for i, struct in enumerate(structures):
            struct.POSCAR.write_frac(dirnames[i] + '/POSCAR')
        self.createSimLauncher(dirnames, launch_command)
        print('Study launcher ready! Hit the magic button:\n\t$ bash launcher.sh')

    @calculate_time
    def createSimLauncher(self, directories, launch_command):
        s = ''
        s += '#!/bin/bash\n'
        s += '\n'
        s += 'for dir in '
        for directory in directories:
            s += '%10s ' % (directory)
        s += '\n'
        s += 'do\n'
        s += '\tcd ${dir}\n'
        s += '\t%s\n' % (launch_command)
        s += '\tcd ../\n'
        s += 'done'

        fopen = open('launcher.sh', 'w')
        fopen.write(s)
        fopen.close()

    @calculate_time
    def study_global_processing(self, infiles, supertag='kbar', subtag='run'):
        # 'infiles' is a list of the filenames to be processed
        namelist, datalist, rootlist = [], [], []
        for i, info in enumerate(infiles):
            names, datasets, rootdirs = self.read_file_study(info, supertag, subtag)
            namelist.append(names)
            datalist.append(datasets)
            rootlist.append(rootdirs)
            sys.stdout.write("Reading files... (%d/%d) \r" % (i+1, len(infiles)))
            sys.stdout.flush()

        sys.stdout.write("\n")
        for i in range(len(infiles)):
           for id_dir, dirs in enumerate(rootlist[i]):
                # If data contains only one row
                if len(datalist[i][id_dir][0,:,0]) == 1:
                    # Initialize this kind of output with the label of the folder as X
                    dataprint = [namelist[i][id_dir]]
                    # Print sequentially all columns of the single row
                    for j in range(len(datalist[i][id_dir][0,0,:])):
                        dataprint.append(datalist[i][id_dir][:, 0, j])
                    header = '# X Y\n'
                else:
                    # If data is a set of XY values, we go into XY1Y2Y3...
                    # and use the names of the subfolders as the info line
                    names = namelist[i][id_dir]
                    s = '# '
                    for ni in range(len(names)):
                        s += '%.2f\t' % (names[ni])
                    header = s + '\n'
                    dataprint = [datalist[i][id_dir][0,:,0]]
                    
                    for j in range(len(namelist[i][id_dir])):
                        dataprint.append(datalist[i][id_dir][j,:,1])
                self.print_to_file(dirs+'/global_' + infiles[i], dataprint, header=header)

        sys.stdout.write("Merging files... %d/%d\r" % (i+1, len(infiles)))
        sys.stdout.flush()
    sys.stdout.write("\n")
