from __future__ import division
from qm.model.model import Model
import numpy as np

class SEA21(Model):
    """ Class for model in SEA21 model BO calculation

        :param object molecule: molecule object
        :param double K: parameter for SEA21 model
        :param double R1: parameter fro SEA21 model
        :param double R2: parameter for SEA21 model
        :param double R3: parameter for SEA21 model
        :param A: parameter for SEA21 model
        :param D: parameter for SEA21  model
        :param G: parameter for SEA21 model
    """
    def __init__(self, molecule, K=0.02, R1=6.0, R2=2.0, R3=3.875, A=3.0, D=0.01, G=0.01):
        # Initialize model common variables
        super(SEA21, self).__init__(None)

        # Define parameters
        self.K = K
        self.R1 = R1
        self.R2 = R2
        self.R3 = R3
        self.A = A
        self.D = D
        self.G = G

        # Set 'l_nacme' with respect to the computational method
        # SEA21 model can produce NACs, so we do not need to get NACME
        molecule.l_nacme = False

        # SEA21 model can compute the gradient of several states simultaneously
        self.re_calc = False

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from SEA21 model BO calculation

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # Initialize diabatic Hamiltonian
        H = np.zeros((2, 2))
        dH = np.zeros((2, 2))
        unitary = np.zeros((2, 2))

        x = molecule.pos[0]

        # Define Hamiltonian
        H[0, 0] = 1. / 2. * self.K * (x - self.R1) ** 2  
        H[1, 1] = 1. / 2. * self.K * (x - self.R2) ** 2 + self.D
        H[1, 0] = self.G * np.exp( -self.A * (x - self.R3) ** 2
        H[0, 1] = H[1, 0]

        # Define a derivative of Hamiltonian
        dH[0, 0] = self.K * (x - self.R1)
        dH[1, 1] = self.K * (x - self.R2)
        dH[1, 0] = - 2. * self.A * self.G * (x - self.R3) * np.exp(- self.A * (x - self.R3) ** 2)
        dH[0, 1] = dH[1, 0]

        # Diagonalization
        a = 4. * H[1, 0] * H[0, 1] + (H[1, 1] - H[0, 0]) ** 2
        sqa = np.sqrt(a)
        tantheta = (H[1, 1] - H[0, 0] - sqa) / H[1, 0]  * 0.5
        theta = np.arctan(tantheta)

        unitary[0, 0] = np.cos(theta)
        unitary[1, 0] = np.sin(theta)
        unitary[0, 1] = - np.sin(theta)
        unitary[1, 1] = np.cos(theta)

        # Extract adiabatic quantities
        molecule.states[0].energy = 0.5 * (H[0, 0] + H[1, 1]) - 0.5 * sqa
        molecule.states[1].energy = 0.5 * (H[0, 0] + H[1, 1]) + 0.5 * sqa

        molecule.states[0].force = np.dot(unitary[:, 1], np.matmul(dH, unitary[:, 1]))
        molecule.states[1].force = np.dot(unitary[:, 0], np.matmul(dH, unitary[:, 0]))

        molecule.nac[0, 1, 0, 0] = np.dot(unitary[:, 0], np.matmul(dH, unitary[:, 1])) / sqa
        molecule.nac[1, 0, 0, 0] = - molecule.nac[0, 1, 0, 0]

