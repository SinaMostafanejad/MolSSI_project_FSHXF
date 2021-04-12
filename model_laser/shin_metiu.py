from __future__ import division
from bo.model.model import Model
import os, shutil, re
import numpy as np

class Shin_Metiu(Model):
    """ Class for Shin-Metiu model BO calculation

        :param object molecule: molecule object
        :param string qm_path: path for Shin-Metiu 1D charge transfer model BO calculation program
    """
    def __init__(self, molecule, qm_path="./"):
        # Initialize model common variables
        super(Shin_Metiu, self).__init__(qm_path)

        # Set 'l_nacme' with respect to the computational method
        # Shin-Metiu model can produce NACs, so we do not need to get NACME
        molecule.l_nacme = False

        # Shin-Metiu model can compute the gradient of several states simultaneously
        self.re_calc = False

    def get_bo(self, molecule, base_dir, istep, bo_list, dt, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from Shin-Metiu BO calculation

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        super().get_bo(base_dir, calc_force_only)
        #self.get_input(molecule,istep)
        self.get_input(molecule,istep,dt)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def get_input(self, molecule,istep,dt):
        """ Generate input file for BO calculation: metiu.in 

            :param object molecule: molecule object
        """
        # Make 'metiu.in' file
        input_shin_metiu = ""
        input_shin_metiu += f"metiu\n{molecule.pos[0, 0]:15.8f}\n{molecule.nst}\n"
        #input_shin_metiu += f"{istep + 1}"
        input_shin_metiu += f"{istep + 1}\n"
        input_shin_metiu += f"{dt}\n"

        # Write 'metiu.in' file
        file_name = "metiu.in"
        with open(file_name, "w") as f:
            f.write(input_shin_metiu)

    def run_QM(self, base_dir, istep, bo_list):
        """ Run Shin-Metiu BO calculation and save the output files to QMlog directory

            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
        """
        # Run Shin-Metiu model
        #qm_command = os.path.join(self.qm_path, "metiu.x")
        qm_command = os.path.join("srun ", self.qm_path, "metiu.x")
        command = f"{qm_command} < metiu.in > metiu.log"
        os.system(command)
        # Copy the output file to 'QMlog' directory
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("metiu.log", os.path.join(tmp_dir, log_step))

    def extract_BO(self, molecule, bo_list, calc_force_only):
        """ Read the output files to get BO information

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # Energy
        file_name = "ENERGY.DAT"
        with open(file_name, "r") as f:
            log_out = f.read()

        if (not calc_force_only):
            for states in molecule.states:
                states.energy = 0.

            tmp_e = '([-]*\S+)'
            energy = re.findall(tmp_e, log_out)
            energy = np.array(energy)
            energy = energy.astype(float)
            for ist in range(molecule.nst):
                molecule.states[ist].energy = energy[ist]

        # Force
        file_name = "FORCE.DAT"
        with open(file_name, "r") as f:
            log_out = f.read()

        if (not calc_force_only):
            for states in molecule.states:
                states.force = np.zeros((molecule.nat, molecule.nsp))

        for ist in bo_list:
            tmp_f = f'{ist + 1:d}\n\s*([-]*\S+E[-]*\S+)'
            force = re.findall(tmp_f, log_out)
            force = np.array(force[0])
            force = force.astype(float)
            force = force.reshape(molecule.nat, molecule.nsp, order='C')
            molecule.states[ist].force = np.copy(force)
		
        # NAC
        file_name = "NAC.DAT"
        with open(file_name, "r") as f:
            log_out = f.read()

        if (not calc_force_only and self.calc_coupling):
            for ist in range(molecule.nst):
                for jst in range(molecule.nst):
                    if (ist == jst):
                        molecule.nac[ist, jst, :, :] = 0.
                    elif (ist < jst):
                        tmp_c = f'{ist + 1:d} {jst + 1:d}\n\s*([-]*\S+)'
                        nac = re.findall(tmp_c, log_out)
                        nac = np.array(nac[0])
                        nac = nac.astype(float)
                        nac = nac.reshape(molecule.nat, molecule.nsp, order='C')
                        molecule.nac[ist, jst] = np.copy(nac)
                    else:
                        molecule.nac[ist, jst] = - molecule.nac[jst, ist]


        # NAC
        file_name = "LASER_COUP.dat"
        with open(file_name, "r") as f:
            log_out = f.read()

        if (not calc_force_only and self.calc_coupling):
            for ist in range(molecule.nst):
                for jst in range(molecule.nst):
                    if (ist == jst):
                        molecule.lasercoup[ist, jst, :, :] = 0.
                    elif (ist < jst):
                        tmp_l = f'{ist + 1:d} {jst + 1:d}\n\s*([-]*\S+)'
                        lasercoup = re.findall(tmp_l, log_out)
                        lasercoup = np.array(lasercoup[0])
                        lasercoup = lasercoup.astype(float)
                        lasercoup = lasercoup.reshape(molecule.nat, molecule.nsp, order='C')
                        molecule.lasercoup[ist, jst] = np.copy(lasercoup)
                    else:
                        molecule.lasercoup[ist, jst] = - molecule.lasercoup[jst, ist]
