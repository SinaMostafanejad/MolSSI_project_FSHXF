rom __future__ import division
from qm.model.model import Model
import numpy as np
from math import erf
from misc import eps, au_to_fs, fs_to_au

class Soft_Core(Model):
    """ Class for 1D soft-core model BO calculation in a real-space grid

        :param object molecule: molecule object
        :param integer nx: the number of grid points
        :param double xmin: lower bound of the 1D space
        :param double xmax: upper bound of the 1D space
        :param double L: the distance between two fixed nuclei
        :param double Rc: the parameter of a moving nucleus
        :param double Rl: the parameter of a fixed nucleus in the left side
        :param double Rr: the parameter of a fixed nucleus in the right side
    """
    def __init__(self, molecule, nx=401, xmin=-20.0, xmax=20.0, E0=0.0093, omega=0.057, tau = 4.8, delta_t=7):
        # Initialize model common variables
        super(Soft_Core, self).__init__(None)

        # Set the grid
        self.nx = nx
        self.xmin = xmin
        self.xmax = xmax


        self.dx = (self.xmax - self.xmin) / np.float(self.nx - 1)
        self.H = np.zeros((self.nx, self.nx))

        #Set the laser parameters
        self.E0 = E0
        self.omega = omega
        self.tau = tau
        self.delta_t = delta_t

        # Set 'l_nacme' with respect to the computational method
        # Soft-Core model can produce NACs, so we do not need to get NACME
        molecule.l_nacme = False

        # Soft-Core model can compute the gradient of several states simultaneously
        self.re_calc = False

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from Soft-Core BO calculation

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # Initialize Hamiltonian
        self.H = 0.
        self.H_laser = 0.

        # Add the kinetic-energy contribution (tridiagonal)
        self.H += - 0.5 * (np.diag([1.] * (self.nx - 1), - 1) + np.diag([- 2.] * self.nx, 0) + \
            np.diag([1.] * (self.nx - 1), 1)) / self.dx ** 2

        self.H_laser += - 0.5 * (np.diag([1.] * (self.nx - 1), - 1) + np.diag([- 2.] * self.nx, 0) + \
            np.diag([1.] * (self.nx - 1), 1)) / self.dx ** 2


 
        x = molecule.pos[0, 0]

        # Add the potential contribution (diagonal)
        xes = [self.xmin + ix * self.dx for ix in range(self.nx)]
        Vs = [self.get_V(x, xe) for xe in xes]
        Vs_laser = [self.get_Vlas(x, xe, istep) for xe in xes]
        
        self.H += np.diag(Vs)
        self.H_laser +=np.diag(Vs_laser)
        
        # Diagonalization
        ws, unitary = np.linalg.eig(self.H)
        ws_laser, unitary = np.lialg.eig(self.H_laser)

        # Sorting eigenvalues in the ascending order and the corresponding eigenvectors
        idx = np.argsort(ws)
        ws = ws[idx]
        unitary = unitary[:, idx]
   
        idx_las = np.argsort(ws_laser)
        ws_laser = ws_laser[idx_las]
        unitary_laser = unitary[:, idx_las]

# Slicing eigenvalues and eigenvectors up to the given number of states
        ws = ws[0:molecule.nst]
        unitary = unitary[:, 0:molecule.nst]
        
        ws_laser = ws_laser[0:molecule.nst]
        unitary_laser = unitary_laser[:, 0:molecule.nst]

        for ist in range(molecule.nst):
            molecule.states[ist].energy = ws[ist]
            molecule.states[ist].energy_laser = ws_laser[ist]

        # Extract adiabatic quantities (BO and QS)
        dVs = [self.get_dV(x, xe) for xe in xes]
        dVijs = np.dot(np.transpose(unitary), np.dot(np.diag(dVs), unitary))
        
        dVs_laser = [self.get_dVlas(x, xe, istep) for xe in xes]
        dVijs_laser = np.dot(np.traspose(unitary), np.dot(np.diag(dVs_laser),uniatry))

        Fs = - np.diag(dVijs)

        Fs_laser = -np.diag(dVijs_laser)

        for ist in range(molecule.nst):
            molecule.states[ist].force = Fs[ist]
            molecule.states[ist].force_laser = Fs_laser[ist]


        for ist in range(molecule.nst):
            for jst in range(ist + 1, molecule.nst):
                molecule.nac[ist, jst, 0, 0] = dVijs[ist, jst] / (ws[jst] - ws[ist])
                molecule.nac[jst, ist, 0, 0] = - molecule.nac[ist, jst, 0, 0]
                molecule.nac_laser[ist,jst, 0, 0] = dVijs_laser[ist, jst] / (ws_laser[jst] - ws_laser[jst])
                molecule.nac_laser[jst, ist, 0, 0] = - molecule.nac_laser[ist, jst, 0, 0]

    def get_V(self, x, xe):
        """ Calculate potential elements of the BO Hamiltonian

            :param double x: the nuclear position
            :param double xe: the electronic position
        """


        V  =  1. / (x + 0.03)**2. - 1. / np.sqrt((xe - x/2.)**2. + 1. ) - \
                1. / np.sqrt((xe + x/2.)**2. + 1. ) 
 

        return V

    def get_dV(self, x, xe):
        """ Calculate del potential elements of the BO Hamiltonian

            :param double x: the nuclear position
            :param double xe: the electronic position
        """
        dV = 0.5 * (x/2. - xe) * ((xe - x /2.)**2. + 1. )**(-1.5) + \
            0.5 * (x/2. + xe) * ((xe + x/2.)**2.+1)**(-1.5) -  \
            x * (x**2 + 0.03)**(-1.5)


        return dV


    def get_Vlas(self, x, xe, istep):
        """ Calculate potential elements of the BO Hamiltonian
 
            :param double x: the nuclear position
            :param double xe: the electronic position
        """
 
 
        Vlas=   1. / (x + 0.03)**2. - 1. / np.sqrt((xe - x/2.)**2. + 1. ) - \
                 1. / np.sqrt((xe + x/2.)**2. + 1. ) + \
                 xe * E0 * exp ((-istep*dt/tau)**2)*cos(omega*istep-deltat)
 
        return Vlas


    def get_dVlas(self, x, xe, istep):
        """ Calculate del potential elements of the BO Hamiltonian
        
            :param double x: the nuclear position
            :param double xe: the electronic position
        """
        dVlas = 0.5 * (x/2. - xe) * ((xe - x /2.)**2. + 1. )**(-1.5) + \
            0.5 * (x/2. + xe) * ((xe + x/2.)**2.+1)**(-1.5) -  \
            x * (x**2 + 0.03)**(-1.5)
            
            
        return dVlas

