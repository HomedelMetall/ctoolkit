from toolkit.global_vars.ext_libs import *

class physicalProperties:
    def __init__(self):
        pass

    def angular_velocities_MD_eig(self, n_vecs, dt):
        # We need error handling in this function
        old_err = np.seterr(all='raise')

        w = np.zeros([len(n_vecs), len(n_vecs[0]), 3], dtype=float)
        print("Processing angular velocities...")

        for step in range(0, len(n_vecs)-1):
            for i in range(len(n_vecs[0])):
                # u, v = n_vecs[step, i], n_vecs[step+1,i]
                u, v = n_vecs[step, i], n_vecs[step+1, i]

                # We consider here the usual derivative averaging
                # in this new context
                w[step, i] = np.cross(u, v)

                # It needs to be normalized
                # But if it's zero, it's zero...
                norm = np.sqrt(np.dot(w[step, i], w[step, i]))
                if not norm < 1E-7:
                    w[step, i] /= np.sqrt(np.dot(w[step, i], w[step, i]))

                # Now we compute w. To do this, we compute the angle
                # between the vectors and how it changes, as if
                # we were sitting in a 2D plane.

                #u, v = n_vecs[step-1, i], n_vecs[step+1,i]

                scalar = np.dot(v,u)
                if scalar > 1.0:
                    theta=np.arccos(1.0)
                else:
                    theta = np.arccos( np.dot(u, v) )#/ (np.sqrt(np.dot(v, v))*np.sqrt(np.dot(u,u))) )

                # We calculate the derivative of the angle to obtain the norm of the angular velocity
                # and apply it to the w vector, taking into account we measure all by
                # double timesteps. Note also that theta= delta theta.
                w[step, i] *= theta/(dt)
            sys.stdout.write("%3.2f...\r" % (100.*(float(step+1))/(len(n_vecs))))
            sys.stdout.flush()

        # We finally treat the boundary cases
        for i in range(len(n_vecs[0])):
            u, v = n_vecs[-1, i], n_vecs[-2, i]
            w[-1, i] = np.cross(u, v)
            norm = np.sqrt(np.dot(w[-1, i], w[-1, i]))
            if not norm < 1E-7:
                w[-1, i] /= np.sqrt(np.dot(w[-1, i], w[-1, i]))

            scalar = np.dot(u,v)
            if scalar > 1.0:
                theta = np.arccos( 1.0 )# / (np.sqrt(np.dot(v, v))*np.sqrt(np.dot(u, u))) )
            else:
                theta = np.arccos( np.dot(u,v) )

            w[-1, i] *= -theta/(dt)

        np.seterr(**old_err)

        return w

    # A function that computes the derivative per step
    # As is, it computes the rotation angle as the crossproduct
    # of two consecutive velocity vectors.
    def angular_velocities_MD_test(self, n_vecs, dt):
        # We need error handling in this function
        old_err = np.seterr(all='raise')

        w = np.zeros([len(n_vecs), len(n_vecs[0]), 3], dtype=float)
        print("Processing angular velocities...")

        for step in range(1, len(n_vecs)-1):
            for i in range(len(n_vecs[0])):
                # u, v = n_vecs[step, i], n_vecs[step+1,i]
                u, v = n_vecs[step, i]-n_vecs[step-1, i], n_vecs[step+1, i]-n_vecs[step,i]

                # We consider here the usual derivative averaging
                # in this new context
                w[step, i] = np.cross(u, v)

                # It needs to be normalized
                # But if it's zero, it's zero...
                norm = np.sqrt(np.dot(w[step, i], w[step, i]))
                if not norm < 1E-7:
                    w[step, i] /= np.sqrt(np.dot(w[step, i], w[step, i]))

                # Now we compute w. To do this, we compute the angle
                # between the vectors and how it changes, as if
                # we were sitting in a 2D plane.
                
                u, v = n_vecs[step-1, i], n_vecs[step+1,i]

                scalar = np.dot(v,u)
                if scalar > 1.0:
                    theta=np.arccos(1.0)
                else:
                    theta = np.arccos( np.dot(u, v) )#/ (np.sqrt(np.dot(v, v))*np.sqrt(np.dot(u,u))) )

                # We calculate the derivative of the angle to obtain the norm of the angular velocity
                # and apply it to the w vector, taking into account we measure all by
                # double timesteps. Note also that theta= delta theta.
                w[step, i] *= theta/(2*dt)
            sys.stdout.write("%3.2f...\r" % (100.*(float(step+1))/(len(n_vecs))))
            sys.stdout.flush()

        """
        # We finally treat the boundary cases
        for i in range(len(n_vecs[0])):
            u, v = n_vecs[-1, i], n_vecs[-2, i]
            w[-1, i] = np.cross(u, v)
            norm = np.sqrt(np.dot(w[-1, i], w[-1, i]))
            if not norm < 1E-7:
                w[-1, i] /= np.sqrt(np.dot(w[-1, i], w[-1, i]))

            scalar = np.dot(u,v)
            if scalar > 1.0:
                theta = np.arccos( 1.0 )# / (np.sqrt(np.dot(v, v))*np.sqrt(np.dot(u, u))) )
            else:
                theta = np.arccos( np.dot(u,v) )

            w[-1, i] *= -theta/(dt)
        """

        np.seterr(**old_err)

        return w[1:-1]

    @calculate_time
    def angular_velocities_MD_GPT(self, n_vecs, dt):
        w = np.zeros([len(n_vecs), len(n_vecs[0]), 3], dtype=float)
        for step in range(1, len(n_vecs)-1):
            u, v = n_vecs[step-1], n_vecs[step+1]
            # we use the numpy cross function to compute cross product 
            # and divide by 2*dt to get the angular velocity
            w[step] = np.cross(u, v) / (2*dt)
        w[0] = w[-1] = w[1] - w[-2]
        return w

    # Calculate angle velocities from the angle description of the system.
    @calculate_time
    def angular_velocities_MD(self, n_vecs, dt):
        # We need error handling in this function
        old_err = np.seterr(all='raise')

        # First we compute the vectorial character of the angular velocity
        # That is, we need to find the vector describing the axis of the rotation.
        # Problem is all this is going to be super slow so I'm gonna try to add some
        # descriptors for the time being.
        w = np.zeros([len(n_vecs), len(n_vecs[0]), 3], dtype=float)
        print("Processing angular velocities...")
        for step in range(1, len(n_vecs)-1):
            for i in range(len(n_vecs[0])):
                u, v = n_vecs[step-1, i], n_vecs[step+1,i]
                # We consider here the usual derivative averaging
                # in this new context
                w[step, i] = np.cross(u, v)
                
                # It needs to be normalized
                try:
                    w[step, i] /= np.sqrt(np.dot(w[step, i], w[step, i]))
                except:
                    pass

                # Now we compute w. To do this, we compute the angle
                # between the vectors and how it changes, as if 
                # we were sitting in a 2D plane.
                try:
                    theta = np.arccos( np.dot(v, u) / (np.sqrt(np.dot(v, v))*np.sqrt(np.dot(u,u))) )
                except:
                    theta = np.arccos(1.0)

                # We calculate the derivative of the angle to obtain the norm of the angular velocity
                # and apply it to the w vector, taking into account we measure all by
                # double timesteps. Note also that theta= delta theta.
                w[step, i] *= theta/(2*dt)
            sys.stdout.write("%3.2f...\r" % (100.*(float(step+1))/(len(n_vecs))))
            sys.stdout.flush()

        # We finally treat the boundary cases
        for i in range(len(n_vecs[0])):
            u, v = n_vecs[1, i], n_vecs[0, i]
            w[0, i] = np.cross(u, v)

            try:
                w[0, i] /= np.sqrt(np.dot(w[0, i], w[0, i]))
            except:
                pass

            try:
                theta = np.arccos( np.dot(v, u) / (np.sqrt(np.dot(v, v))*np.sqrt(np.dot(u,u))) )
            except:
                theta = np.arccos(1.0)
            w[0, i] *= -theta/(dt)
        
        for i in range(len(n_vecs[0])):
            u, v = n_vecs[-1, i], n_vecs[-2, i]
            w[-1, i] = np.cross(u, v)
            try:
                w[-1, i] /= np.sqrt(np.dot(w[-1, i], w[-1, i]))
            except:
                pass

            try:
                theta = np.arccos( np.dot(v, u) / (np.sqrt(np.dot(v, v))*np.sqrt(np.dot(u,u))) )
            except:
                theta = np.arccos(1.0)
            w[-1, i] *= -theta/(dt)

        np.seterr(**old_err)

        return w


    def angular_velocities_MD_deprecated(self, angles, dt):
        angular_velocities = np.copy(angles)
        angular_velocities[0] = (angles[1] - angles[0])/dt
        angular_velocities[-1] = (angles[-1] - angles[-2])/dt
   
        post = np.copy(angles[2:])
        pre = np.copy(angles[:-2])
        angular_velocities[1:-1] = (post-pre)/(2*dt)

        #for i in range(1,len(angles)-1):
        #    angular_velocities[i] = (angles[i+1] - angles[i-1])/(2*dt)

        return angular_velocities

    @calculate_time
    def calculate_angle_full(self, v1, v2, e_crossproduct):
        cosphi = np.inner(v1, v2) / (np.linalg.norm(v1)*np.linalg.norm(v2))
        sinphi = np.dot(e_crossproduct, np.cross(v1,v2))
        angle = 180*np.arccos(cosphi)/np.pi
        if sinphi>0:
            return angle
        else:
            return 360-angle

    # e_crossproduct is the e vector normal to the plane v
    # We could compute this internally but then we would need
    # to know which vector projects into which vector, v1, into v2.
    @calculate_time
    def calculate_angle_full_original(self, v1, v2, e_crossproduct):
        cosphi = np.dot(v1, v2)
        sinphi = np.dot(e_crossproduct, np.cross(v1, v2))
        acos = 180*np.arccos(cosphi)/np.pi
        if sinphi >= 0.:
            angle = acos
        else: 
            #sinphi < 0.:
            angle = 180+(180-acos)
        return angle

    @calculate_time
    def calculate_angle(self, v1, v2):
        angle = np.arccos(np.dot(v1, v2)/np.sqrt(np.dot(v1, v1)*np.dot(v2, v2))) 
        return angle

    @calculate_time
    def calculate_polarization(self, born_charges, ref_struct, final_struct):
        displacements = np.copy(final_struct.at_cart - self.frac_to_cart(final_struct.B, ref_struct.at_frac)) # Any cell basis! :)

        P = np.zeros([3])
        for alpha in range(3):
            for atom in range(len(born_charges)): # this could be also final_struct.at_frac
                for beta in range(3):
                    P[alpha] += (born_charges[atom, beta, alpha]*displacements[atom, beta]) / final_struct.volume

        # Conversion factors
        EtoC = 1.6021766E-19
        AtoM = 1.0E-20
        
        P *= EtoC*(1.0E6)/(AtoM*10000.0)
        print("Polarization (uC/cm^2): %.8f %.8f %.8f" %(P[0], P[1], P[2]))
        print("Modulus (uC/cm^2): %.8f" % (np.sqrt(np.dot(P, P))))

        return P

    # ATTENTION: we assume SCALEUP formalism, in which cells *are* defined.
    @calculate_time
    def calculate_polarization_direction(self, born_charges, ref_struct, final_struct):
        # Conversion factors
        EtoC = 1.6021766E-19
        AtoM = 1.0E-20

        # Get displacements from RS
        displacements = np.copy(final_struct.at_cart - self.frac_to_cart(final_struct.B, ref_struct.at_frac)) # Any cell basis! :)
        
        # Get cell ordering
        #C[i,j,k]->displacements[atom, beta] <- ordering
        C = np.copy(ref_struct.cell_id)

        # Get polarization of each cell
        nx, ny, nz = len(C), len(C[0]), len(C[0,0]) # Not sure it works?
        Px = np.zeros([nx, 3], dtype=float)
        Py = np.zeros([ny, 3], dtype=float)
        Pz = np.zeros([nz, 3], dtype=float)
        P = np.zeros([len(C), len(C[0]), len(C[0,0]), 3], dtype=float)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    # Redefine displacements into disp <- according to C[i,j,k]
                    for alpha in range(3):
                        for atom in range(len(born_charges)): # here we have to go to 5 atoms cell?
                            for beta in range(3):
                                P[i,j,k,alpha] += (born_charges[atom, beta, alpha]*displacement[atom, beta]) / (final_struct.volume/(nx*ny*nz))

        # Unit conversion
        P *= EtoC*(1.0E6)/(AtoM*10000.0)

        # Get polarization per direction
        for i in range(nx):                       
            for j in range(ny):
                for k in range(nz):
                    Px[i] += P[i,j,k]
                    Py[j] += P[i,j,k]
                    Pz[k] += P[i,j,k]

        self.print_to_file('P_along_x.dat', [Px])
        self.print_to_file('P_along_y.dat', [Py])
        self.print_to_file('P_along_z.dat', [Pz])

    @calculate_time
    def BMprime(self, V, e0,v0,k,ep): 

        return (160.217646)* (1.5* v0 * k)* ( (1.0+ 2.0*ep)*(v0/V)**(4.0/3.0)*(1.0/V) - (ep*(v0/V)**(2.0)*(1.0/V)) - (1.0 + ep)*(v0/V)**(2.0/3.0)*(1.0/V) )

    @calculate_time
    def BM(self, V, e0,v0,k,ep):

        return e0 + (1.5 * v0 * k) * ( 0.75*(1.0 + (2.0*ep)) * (v0/V)**(4.0/3.0) - (0.5*ep)*(v0/V)**2.0 - 1.5*(1.0 + ep)*(v0/V)**(2.0/3.0) + 0.5*(ep + 1.5) )

    @calculate_time
    def GF(self, init_params, VE):
        GFvalue = 0.0
        [e0,v0,k,ep] = init_params
        for i in range(len(VE)):
            GFvalue += (VE[i,1]-self.BM(VE[i,0], e0, v0, k, ep))**2/len(VE)
        return GFvalue

    @calculate_time
    def fit_BM(self, V, E):
        init_params = [-59, 156, 1.34, 1.1]
        VE = np.zeros([len(V),2], dtype=float)
        for i in range(len(V)):
            VE[i,0] = V[i]+0.0
            VE[i,1] = E[i]+0.0

        result1 =  optimize.minimize(self.GF, init_params, args=(VE))['x']
        return optimize.minimize(self.GF, result1, args=(VE))

    @calculate_time
    def PV_range(self, range_V, BMparams, filename=None):
        if filename is None: filename = 'PV.dat'
        fopen = open(filename, 'w')
        s = '#V(Angstrom**3)\tP(GPa)\n'
        numpoints = 100
        for i in range(numpoints+1):
            V = range_V[0] + (float(i)/float(numpoints))*(range_V[1]-range_V[0])
            P = self.BMprime(V, BMparams[0], BMparams[1], BMparams[2], BMparams[3])
            s += '%.8e\t%.8e\n' % (V, P)
        fopen.write(s)
        fopen.close()    

    @calculate_time
    def EV_range(self, range_V, BMparams, filename=None):
        if filename is None: filename = 'EV.dat'
        fopen = open(filename, 'w')
        s = '#V(Angstrom**3)\tE(eV)\n'
        numpoints = 100
        for i in range(numpoints+1):
            V = range_V[0] + (float(i)/float(numpoints))*(range_V[1]-range_V[0])
            E = self.BM(V, BMparams[0], BMparams[1], BMparams[2], BMparams[3])
            s += '%.8e\t%.8e\n' % (V, E)
        fopen.write(s)
        fopen.close()

    # It is prepared to call as process_EV([VASP1.outcar, VASP2.outcar, ...])
    @calculate_time
    def process_EV(self, V, E):
        E_array = np.array(E, dtype=float)
        V_array = np.array(V, dtype=float)

        fit_result = self.fit_BM(V_array, E_array)
        # Process fit
        BMparams = np.copy(fit_result['x'])
        testparams = [-59.6209, 157.55, 0.0829353, 0.344118]
        self.EV_range([np.min(V_array), np.max(V_array)], BMparams, 'EV_fit.dat')
        self.PV_range([np.min(V_array), np.max(V_array)], BMparams, 'PV_fit.dat')

        fopen = open('EV_raw.dat', 'w')
        s = '#V (Angstrom**3)\tE(eV)\n'
        for i in range(len(V)):
            s += '%.8e\t%.8e\n' % (V[i], E[i])
        fopen.write(s)

        return BMparams

    @calculate_time
    def solve_VP(self, final_pressure, BMparams, verbose=False):
        FP = final_pressure
        [e0, v0, k, ep] = BMparams
        tol = 1E-7

        # First sample
        P = 0
        for i in range(1, 2000):
            last_step = FP-P
            P = self.BMprime(float(i), e0, v0, k, ep)+0.0
            this_step = FP - P
        
            if(last_step<0.0 and this_step>0.0):
                init_range = [float(i-1), float(i)]

        Pdiff = 1.0
        [V1, V2] = init_range
        steps = 0
        while(Pdiff>tol):
            new_volume = (V1+V2)/2
            P = self.BMprime(new_volume, e0, v0, k, ep)
            if P < final_pressure:
                [V1, V2] = [V1, new_volume]
            if P > final_pressure:
                [V1, V2] = [new_volume, V2]
            steps += 1
            Pdiff = np.abs(P-FP)

        if verbose: print("Volume(Angstrom**3) and pressure(GPa) (%dsteps):"%(steps), new_volume, P)
        return new_volume

    @calculate_time
    def Fvib_from_PDOS(self, omegaTHZ, fft, T):
        # Frequency in omega is expected in THz
        kB = 8.617333262E-5 #eV*K^-1
        h = 4.135667696E-15 # eV*s
        omega = 1E12*omegaTHZ
        integrand = np.zeros([int(len(omega))])
        integrand_sinh = np.zeros([int(len(omega))])
        for i in range(1, int(len(omega))):
            integrand[i] = fft[i]*(kB*T*np.log(1.0-np.exp(-(h*omega[i])/(kB*T)))+h*omega[i]/2.0) # VERSION         1
            integrand_sinh[i] = kB*T*(fft[i]*np.log(2.*np.sinh(h*omega[i]/(2.*kB*T))))

        #print(omega,fft)
        F = self.func_integral(omega, integrand_sinh)
        #print(F)
        return F

    # Direct extract of G from PDOS
    @calculate_time
    def S_from_PDOS(self, omegaTHZ, fft, T):
        # Frequency in omega is expected in THz
        kB = 8.617333262E-5 #eV*K^-1
        h = 4.135667696E-15 # eV*s
        omega = 1E12*omegaTHZ
        integrand = np.zeros([int(len(omega))])
        for i in range(1, int(len(omega))):
            argument = h*omega[i]/(2.*kB*T)
            
            integrand[i] = fft[i]*(kB*np.log(2*np.sinh(argument)) - (h*omega[i]*(np.cosh(argument)/np.sinh(argument)))/(2*T))

        S = self.func_integral(omega, integrand)

        return S

    @calculate_time
    def S_from_2Dhist(self, xedges, yedges, hist):
        S = .0
        # There is no kbar...
        old_err = np.seterr(all='raise')
        np.seterr(all = 'raise')
        number_of_bins = len(hist)*len(hist[0])
        norm_hist = np.copy(hist)
        xlen = xedges[-1] - xedges[0]
        ylen = yedges[-1] - yedges[0]
        area_element = xlen*ylen/number_of_bins
        for i in range(len(norm_hist)):
            for j in range(len(norm_hist[i])):
                try:
                    S += area_element*norm_hist[i,j]*np.log(area_element*norm_hist[i,j])
                except:
                    continue

        np.seterr(**old_err)
        S = S - np.log(number_of_bins)
        return S

    @calculate_time
    def C_from_PDOS(self, omega, fft, T):
        kB = 8.617333262E-5 #eV*K^-1
        hbar = 6.582119569E-16 # eV*s
        integrand = np.zeros([int(len(omega))])
        for i in range(1, int(len(omega))):
            argument = hbar*1.0E12*omega[i]/(kB*T)

            integrand[i] = fft[i]*kB*(argument**2)*(np.exp(argument)/((np.exp(argument)-1)**2))

        C = self.func_integral(1E12*omega, integrand)

        return C

    @calculate_time
    def G_from_Fvib(self, Fvib, EP, press, vol):
        # Assumes the input units are coherent... :)
        G = EP + Fvib + (press*vol)
        return G

    # Just a wrapper for np.fit
    @calculate_time
    def fit_func(self, f, xdata, ydata):
        popt, pcov = optimize.curve_fit(f, xdata, ydata)
        return popt

    @calculate_time
    def func_integral(self, x, y):
        return np.trapz(y, x)

    @calculate_time
    def func_interpol(self, x, y):
        # We multiply by N the number of points
        n = 10
        new_x = np.zeros([(len(x)-1)*n], dtype=float)
        for i in range(len(x)-1):
            for j in range(n):
                new_x[j+n*i] = x[i] + (x[i+1]-x[i])*(float(j)/float(n))

        return new_x, np.interp(new_x, x, y)

    @calculate_time
    def func_spline(self, x, y):
        # We multiply by 3 the number of points, for example
        n = 5 
        new_x = np.zeros([(len(x)-1)*5], dtype=float)
        for i in range(len(x)-1):
            for j in range(5):
                new_x[j+5*i] = x[i] + (x[i+1]-x[i])*(float(j)/float(n))

        tck = interpolate.splrep(x, y)

        y_spline = interpolate.splev(new_x, tck)

        return new_x, y_spline

    @calculate_time
    def func_derivative(self, x, y):
        # We multiply by 3 the number of points, for example
        n = 3
        new_x = np.zeros([(len(x)-1)*3], dtype=float)
        for i in range(len(x)-1):
            for j in range(3):
                new_x[j+3*i] = x[i] + (x[i+1]-x[i])*(float(j)/float(n))

        tck = interpolate.splrep(x, y)

        y_der = interpolate.splev(new_x, tck, der=1)

        return new_x, y_der

    @calculate_time
    def func_derivative_direct(self, x, y):
        dx = np.zeros([len(x)-2], dtype=float)
        for i in range(1,len(x)-2):
            dx[i] = (y[i+1]-y[i]) / (x[i+1]-x[i])
            dx[i] += (y[i]-y[i-1]) / (x[i]-x[i-1])
            dx[i] *= 1/2.

        return x[1:-1], dx

    # ACFunc is the autocorrelation function
    # dt is the spacing at which <v(0)*v(dt)> is computed, in seconds
    # Result omega is in THz, and vacf is unbiased.
    # Returns omega, fft
    @calculate_time
    def PDOS_from_ACF(self, ACF, dt, numatoms, normalize=True):
        N=len(ACF)
        omegalong = np.fft.fftfreq(2*N-1, dt)*1E-12
        PDOSlong = np.abs(np.fft.fft(ACF-np.average(ACF), 2*N-1))

        # Only return the non-negative-omega spectrum
        omega = np.copy(omegalong[:int(len(omegalong)/2)])
        PDOS = np.copy(PDOSlong[:int(len(omegalong)/2)])

        if normalize==True:
            norm = self.func_integral(1.0E12*omega, PDOS)
            normalized_PDOS = 3*numatoms*PDOS/norm
            result_PDOS = np.copy(normalized_PDOS)
        else:
            result_PDOS = np.copy(PDOS)

        return omega, result_PDOS

    @calculate_time
    def autocorrelate_vector(self, X):
        result = signal.correlate(X, X, mode='full', method='fft')
        return result[int(result.size/2):]

    @calculate_time
    def autocorrelate_peratom(self, window, i, shVacf):
        shVacf[0] += self.autocorrelate_vector(window[:, i, 0])
        shVacf[1] += self.autocorrelate_vector(window[:, i, 1])
        shVacf[2] += self.autocorrelate_vector(window[:, i, 2])

    @calculate_time
    def autocorrelate_window(self, window, shVacf):
        num_atoms = len(window[0])
        for i in range(num_atoms):
            shVacf[0] += self.autocorrelate_vector(window[:, i, 0])
            shVacf[1] += self.autocorrelate_vector(window[:, i, 1])
            shVacf[2] += self.autocorrelate_vector(window[:, i, 2])

    @calculate_time
    def compute_VACF_parallel(self, v_array, num_windows):
        num_steps = len(v_array)
        num_atoms = len(v_array[0])

        # Parallel routine initialization
        manager = multiprocessing.Manager()
        shVacf = manager.list([manager.list([0.0]*int(num_steps/2))]*3)

        nprocs = 16
        parsteps = int(num_windows/nprocs)
        parsteps_left = int(num_windows%nprocs)

        vacf = np.zeros([int(num_steps/2)], dtype=float)
        vacfx = np.zeros([int(num_steps/2)], dtype=float)
        vacfy = np.zeros([int(num_steps/2)], dtype=float)
        vacfz = np.zeros([int(num_steps/2)], dtype=float)
        
        for par in range(parsteps):
            runners = []
            for proc in range(nprocs):
                i_window = int(par*nprocs + proc)
                window = np.copy(v_array[i_window:int(i_window+num_steps/2)])
                Args = (window, shVacf)
                runners.append(multiprocessing.Process(
                        target=self.autocorrelate_window,
                                                  args=Args))

            #Launch processes
            for p in runners:
                p.start()
            for p in runners:
                p.join()
            
        runners = []
        for par in range(parsteps_left):
            i_window = int(par + parsteps*nprocs)
            window = np.copy(v_array[i_window:int(i_window+num_steps/2)])
            Args = (window, shVacf, )
            runners.append(multiprocessing.Process(
                    target=self.autocorrelate_window,
                                                args=Args))
        #Launch processes
        for p in runners:
            p.start()
        for p in runners:
            p.join()

        vacf = np.array(shVacf, dtype=float)
        vacf = np.copy((vacf[0] + vacf[1] + vacf[2])/(3.0*num_windows*num_atoms))

        del shVacf
        manager.shutdown()

        return np.copy(vacf)

    @calculate_time
    def compute_VACF_serial(self, v_array, num_windows):
        num_steps = len(v_array)
        num_atoms = len(v_array[0])

        vacf = np.zeros([int(num_steps/2)], dtype=float)
        vacfx = np.zeros([int(num_steps/2)], dtype=float)
        vacfy = np.zeros([int(num_steps/2)], dtype=float)
        vacfz = np.zeros([int(num_steps/2)], dtype=float)
        for i in range(num_windows):
            sp = int(i*(num_steps/2)/num_windows)
            ep = int(sp + num_steps/2)
            window = np.copy(v_array[sp:ep])

            # Arizona contribution:
            for j in range(num_atoms):
                vacfx += self.autocorrelate_vector(window[:, j, 0])
                vacfy += self.autocorrelate_vector(window[:, j, 1])
                vacfz += self.autocorrelate_vector(window[:, j, 2])

        vacf = np.copy((vacfx + vacfy + vacfz)/(3.0*num_windows*num_atoms))

        return np.copy(vacf)
    
    @calculate_time
    def compute_VACF_direct(self, v_array, num_windows):
        num_steps, num_atoms, _ = v_array.shape
        # compute the VACF for each window
        windows = np.array_split(v_array[:num_steps//2], num_windows)
        vacf = np.sum(np.sum(w[:,:,0]*w[0,:,0]+w[:,:,1]*w[0,:,1]+w[:,:,2]*w[0,:,2], axis=1) for w in windows)
        # normalize the result
        vacf /= (3.0*num_windows*num_atoms)
        return vacf

    # Wrapper for parallelization
    # Inner parallelization: parallellization of the v(0)v(t) product, for several t simultaneously. This is the most efficient.
    # Outer parallelization: parallelize the windows.
    # Serial: no parallellization
    @calculate_time
    def compute_VACF(self,  v_array, num_windows, mode='direct'):

        if len(v_array)/2 < num_windows:
            print("The windowing mode requires that you "+ 
                    "have less windows than half the number steps! Aborting.")
            sys.exit()

        if (mode != 'direct') and (mode != 'serial') and (mode != 'parallel'):
            print("Default VACF calculation mode: serial")
            mode = 'direct'

        if mode == 'parallel':
            result = self.compute_VACF_parallel(v_array, num_windows)
        if mode == 'serial':
            result = self.compute_VACF_serial(v_array, num_windows)
        if mode == 'direct':
            result = self.compute_VACF_direct(v_array, num_windows)

        return result

    # We have to see if np.dot(vec of vecs) is what we want...
    # ... no, it doesn't. Hello.
    # So what we have now is a self-made RACF calculation that should give the right
    # result whatsoever...
    @calculate_time
    def compute_RACF(self, v_array, num_windows=100, normalize=True):
        num_steps = len(v_array)
        num_vectors = len(v_array[0])
        if normalize == True:
            n_array = self.normalize_vectors_MD(v_array)
        else:
            n_array = np.copy(v_array)

        racf = np.zeros([int(num_steps/2)], dtype=float)
        for i in range(num_windows):
            start = int(i*(num_steps/2)/num_windows)
            window = np.copy(n_array[start:int(start+num_steps/2)])
            for j in range(len(window)):
                racf[j] += np.sum(np.einsum('...j,...j', window[0], window[j]))

        return np.copy(racf / (num_windows*num_vectors)) 

    def ACF_2D(self, v_array, num_windows=100):
        num_steps = len(v_array)
        num_vectors = len(v_array[0])

        acf = np.zeros([int(num_steps/2)], dtype=float)
        acf2 = np.copy(acf)
        for i in range(num_windows):
            start = int(i*(num_steps/2)/num_windows)
            window = np.copy(v_array[start:int(start+num_steps/2)])
            for j in range(len(window)):
                acf[j] += np.sum(np.einsum('...j,...j', window[0], window[j]))
                for k in range(len(window[i])):
                    acf2[j] += np.dot(window[0,k], window[j,k])

        self.print_to_file('test_ACF.dat', [acf])
        self.print_to_file('test_ACF.dat', [acf2])
        return np.copy(acf / (num_windows*num_vectors))


    @calculate_time
    def compute_RACF_test(self, v_array, normalize=True):
        num_steps = len(v_array)
        num_vectors = len(v_array[0])
        if normalize == True:
            n_array = self.normalize_vectors_MD(v_array)
        else:
            n_array = np.copy(v_array)

        racf = np.zeros([int(num_steps)], dtype=float)
        norm_ct = np.zeros([int(num_steps)], dtype=float)
        for t1 in range(num_steps):
            for t2 in range(t1, num_steps):
                for k in range(len(n_array[0])):
                    racf[t2-t1] += np.dot(n_array[t1, k], n_array[t2, k])
                    norm_ct[t2-t1] += 1

        for i in range(num_steps):
            racf[i] *= 1./norm_ct[i]

        return np.copy(racf)

    # ... and on the other hand, an Arizona-like function
    # that should also provide a result. To be compared really.
    @calculate_time
    def compute_RACF_serial(self, v_array, num_windows=100, normalize=True):
        num_steps = len(v_array)
        num_vecs = len(v_array[0])
        if normalize == True:
            n_array = self.normalize_vectors_MD(v_array)
        else:
            n_array = np.copy(v_array)

        racf = np.zeros([int(num_steps/2)], dtype=float)
        racfx = np.zeros([int(num_steps/2)], dtype=float)
        racfy = np.zeros([int(num_steps/2)], dtype=float)
        racfz = np.zeros([int(num_steps/2)], dtype=float)
        for i in range(num_windows):
            sp = int(i*(num_steps/2)/num_windows)
            ep = int(sp + num_steps/2)
            window = np.copy(n_array[sp:ep])

            # Arizona contribution:
            for j in range(num_vecs):
                racfx += self.autocorrelate_vector(window[:, j, 0])
                racfy += self.autocorrelate_vector(window[:, j, 1])
                racfz += self.autocorrelate_vector(window[:, j, 2])

        racf = np.copy((racfx + racfy + racfz)/(num_windows*num_vecs))

        return np.copy(racf)

