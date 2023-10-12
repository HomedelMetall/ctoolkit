from toolkit.global_vars.ext_libs import *

class manageDirs:
    def __init__(self):
        self.getDirs()
        self.getCalculationDirs()

    @calculate_time
    def fileExists(self, filename):
        try:
            open(filename)
            return 0
        except:
            return -1

    @calculate_time
    def getDirs(self, path='./'):
        list_of_dirs = sorted(os.listdir(path))

        return list_of_dirs

    # spanDirs respects order of variable! :)
    @calculate_time
    def spanDirs(self, var, header=None, style='int'):
        if header is None: header = ''
        dirnames = []
        for i in range(len(var)):
            if style=='float': dirname = header + '%.4f' % (var[i])
            elif style=='exp': dirname = header + '%.4e' % (var[i])
            elif style=='int': dirname = header + '%d' % (var[i])
            try:
                os.mkdir(dirname)
            except:
                print('WARNING! Directory %s already exists! Careful overwritting...' % (dirname))
                
            dirnames.append(dirname)

        return dirnames

    # This is just for VASP atm
    @calculate_time
    def getCalculationDirs(self):
        list_calc_dirs = []
        for directory in self.getDirs():
            filename = directory+'/OUTCAR'
            if self.fileExists(filename) == 0:
                list_calc_dirs.append(directory)

        return list_calc_dirs

    @calculate_time
    def copyBASEdir(self, directories):
        if os.path.isdir('INPUT_BASE'):
            for directory in directories:
                os.system('cp INPUT_BASE/* %s/' % (directory))
        else:
            print('You need to create a INPUT_BASE in your root folder \nand put there all that you want to be copied!')
            return -1
        return 0

