from molecule import Molecule
import qm, mqc
from misc import data

data[f"X1"] = 1836 # au

geom = """
1
Soft-Core model
X1       2.0     0.0
"""

mol = Molecule(geometry=geom, nstates=2, ndim=1, unit_pos='au', l_model=True)

qm = qm.model.Soft_Core(molecule=mol)

md = mqc.SH(molecule=mol , nsteps=3500, dt=0.01, unit_dt='fs', istate=1)

#md = mqc.SHXF(molecule=mol, nsteps=2890, nesteps=1, dt=0.5, unit_dt='au', wsigma=0.1, istate=1, propagation="density")

md.run(qm=qm, output_dir=f"./TRAJ_SC_test", l_save_scr=True, l_save_qm_log=False)
