from ctoolkit.global_vars.ext_libs import *
import multiprocessing as mp

# Timers may not be properly displayed for everything that happens
# behind the scenes of the parallel processing here. For instance,
# init_pool_processes, which is the initialization of the process
# inside its own codespace, will not be accounted for the timers
# decorator. This could be fixed using shared memory variables but
# it seems a bit too much just to measure the cost of initializing
# some data locks...
class parallelSuite:
    @calculate_time
    def parallel_manager(self):
        self.pmanager = mp.Manager()
        self.pool = self.create_pool(self.get_lock())
        return self.pmanager

    def get_lock(self):
        self.mpLock = mp.Lock()
        return self.mpLock
    
    def mplist(self, shape=None):
        if self.pmanager is None:
            print("Can't create a list without initializing a manager!")
            return None
        if shape is None:
            return self.pmanager.list()
        else:
            return self.pmanager.list(shape)
    
    @calculate_time
    def create_pool(self, data_lock):
        return mp.Pool(initializer=self.init_pool_processes,
                       initargs=(data_lock, ))

    @calculate_time
    def init_pool_processes(self, the_lock):
        # Init the global locks
        global mpLock
        mpLock = the_lock

    @calculate_time
    def close_pool(self, pool):
        pool.close()

    #@calculate_time
    def run_parallel(self, func, iter_args, chksz = None, pool = None):
        # Check that pool has been initiated.
        # Else, instantiate a new pool
        if self.pool is None:
            self.parallel_manager()
            print("Initialized PM")
        # If pool != None
        if not pool is None:
            self.close_pool(self.pool)
            self.pool = pool
            print("External pool used")

        # Look for chunksize definition
        if chksz is None:
            chksz = 1

        # A bit of an idea on how to code for this.
        # One might be in one of two scenarios.
        # Scenario 1: we use the result of the function.
        # ===========
        #  
        # When parallelizing, "func" here should be
        # representing the step of an iteration.
        # "func" then can be built to return a certain
        # value or object, which would be returned from here.
        # Then the user can use the result computed in parallel
        # to store however the data.
        # 
        # Scenario 2: we manage an external array in parallel
        # ===========
        #
        # Imagine that each step consists of adding an item
        # to a list array. This, which would happen in a
        # *DISORDERED* way (async is applied), could be coded
        # by creating an external variable with parallelSuite.pmlist()
        # and passing it to the the function in iter_args.
        # Additionally, one can even obtain an ordered result
        # if additional arguments are considered, as the
        # index of the iterator being provided as well in iter_args
        # for the function to take them into account for given-size
        # arrays (parallelSuite.pmlist(shape=[[]*N]*M))
        # In that case, you will also need to provide the data
        # lock and the pool externally:
        #
        #   lock = tools.get_lock()
        #   pool = tools.create_pool(lock)
        #   result = tools.run_parallel(func, iterargs, pool=pool)
        #
        # Your "func" then should write the data into the array
        # in a way following this idea:
        #
        #   with lock:
        #       shared_array.append(whatever_data)


        result = self.pool.starmap_async(func, iter_args, chunksize=chksz)
        result.wait()
        
        # Shape of AsyncResult Objects is ugly.
        # So we turn it into a nicer, more handable list object
        nice_result = []
        for value in result.get():
            nice_result.append(value)

        return nice_result

    

