from __future__ import division
from qm.model.model import Model
import numpy as np
from math import erf
from misc import eps, au_to_fs, fs_to_au, nm_to_au

class Soft_Core(Model):
    """ Class for 1D soft-core model BO calculation in a real-space grid

        :param object molecule: molecule object
        :param integer nx: the number of grid points
        :param double xmin: lower bound of the 1D space
        :param double xmax: upper bound of the 1D space
        :param alpha: laser on/off
        :param double E0: laser amplitude
        :param double omega: frequency in au
        :param double tau: pulse duration
        :param double delta_t: ttime delay
    """
    def __init__(self, molecule, unit_dt, nx=401, xmin=-20.0, xmax=20.0, E0=0.0093, omega=800, tau = 4.8, delta_t=7.0, alpha=0.0):
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
        self.omega = nm_to_au / omega  
        self.alpha = alpha
     
        ##Conversion unit laser parameters
        if (unit_dt == 'au'):
            self.delta_t = delta_t
            self.tau = tau
        elif (unit_dt == 'fs'):
            self.delta_t = delta_t * fs_to_au
            self.tau = tau * fs_to_au


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

        # Add the kinetic-energy contribution (tridiagonal)
        self.H += - 0.5 * (np.diag([1.] * (self.nx - 1), - 1) + np.diag([- 2.] * self.nx, 0) + \
            np.diag([1.] * (self.nx - 1), 1)) / self.dx ** 2

 
        x = molecule.pos[0, 0]

        # Add the potential contribution (diagonal)
        xes = [self.xmin + ix * self.dx for ix in range(self.nx)]
        Vs = [self.get_V(x, xe, istep, dt) for xe in xes]
        
        self.H += np.diag(Vs)
        
        # Diagonalization
        ws, unitary = np.linalg.eig(self.H)

        # Sorting eigenvalues in the ascending order and the corresponding eigenvectors
        idx = np.argsort(ws)
        ws = ws[idx]
        unitary = unitary[:, idx]

# Slicing eigenvalues and eigenvectors up to the given number of states
        ws = ws[0:molecule.nst]
        unitary = unitary[:, 0:molecule.nst]

        for ist in range(molecule.nst):
            molecule.states[ist].energy = ws[ist]

        # Extract adiabatic quantities (BO and QS)
        dVs = [self.get_dV(x, xe) for xe in xes]
        dVijs = np.dot(np.transpose(unitary), np.dot(np.diag(dVs), unitary))
#        dlaser_s = [self.get_dlaser(x, xe, istep, dt) for xe in xes]
#        dVijs_laser = np.dot(np.transpose(unitary), np.dot(np.diag(dlaser_s), unitary))
        
        Fs = - np.diag(dVijs)

        for ist in range(molecule.nst):
            molecule.states[ist].force = Fs[ist]

        for ist in range(molecule.nst):
            for jst in range(ist + 1, molecule.nst):
                molecule.nac[ist, jst, 0, 0] = dVijs[ist, jst] / (ws[jst] - ws[ist])
                molecule.nac[jst, ist, 0, 0] = - molecule.nac[ist, jst, 0, 0]
#                molecule.laser_coup[ist, jst, 0, 0] = dVijs_laser [ist, jst] / (ws[jst] - ws[ist]) 
#                molecule.laser_coup[jst, ist, 0, 0] = - molecule.laser_coup[ist, jst, 0, 0]

    def get_V(self, x, xe, istep, dt):
        """ Calculate potential elements of the BO Hamiltonian

            :param double x: the nuclear position
            :param double xe: the electronic position
            :param alpha: laser on/off
            :param istep: time step
        """


        V  =  1. / (x + 0.03)**2. - 1. / np.sqrt((xe - x/2.)**2. + 1. ) - \
                1. / np.sqrt((xe + x/2.)**2. + 1. ) + \
            + self.alpha *  xe * self.E0 * np.exp ((-istep*dt-self.delta_t/self.tau)**2)*np.cos(self.omega*istep-self.delta_t)  

        return V

    def get_dV(self, x, xe):
        """ Calculate potential elements of the BO Hamiltonian

            :param double x: the nuclear position
            :param double xe: the electronic position
        """
        dV = 0.5 * (x/2. - xe) * ((xe - x /2.)**2. + 1. )**(-1.5) + \
            0.5 * (x/2. + xe) * ((xe + x/2.)**2.+1)**(-1.5) -  \
            x * (x**2 + 0.03)**(-1.5)

        return dV

    def get_dlaser(self, x, xe, istep, dt):
        """ Calculate laser derivative

          :param double x:the nuclear position
          :param double xe: the electronic position
          :param istep: time step
          :param E0: laser amplitude
          :param omega: laser frequency
          :param tau: pulse width
          :param delta_t: time delay
          :param alpha: laser on/off

        """
        dlaser = self.alpha * (-xe * self.E0 * self.omega * np.exp ((-istep*dt-self.delta_t/self.tau)**2.) * np.sin(self.omega*(istep*dt-self.delta_t)) - \
                2. * self.E0 * xe * (istep*dt - self.delta_t) * np.exp ((-istep*dt-self.delta_t/self.tau)**2.) * np.cos (self.omega * (istep*dt-self.delta_t)) / self.tau**2.)

        return dlaser
