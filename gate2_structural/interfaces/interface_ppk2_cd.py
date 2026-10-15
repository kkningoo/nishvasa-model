from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.core.kinematics import MoveMap

init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false")
sf = get_score_function()

def analyse(pdb, label):
    pose = pose_from_pdb(pdb)
    print(f"\n=== {label} ===")
    print(f"  score before relax {sf(pose):10.1f}")
    mm = MoveMap(); mm.set_bb(True); mm.set_chi(True)
    fr = FastRelax(sf, 5); fr.set_movemap(mm); fr.apply(pose)
    print(f"  score after relax  {sf(pose):10.1f}")
    ia = InterfaceAnalyzerMover(1, False, sf, True, False, True, False)
    ia.apply(pose)
    d = ia.get_all_data()
    dg, sasa = ia.get_interface_dG(), ia.get_interface_delta_sasa()
    print(f"  dG_separated       {dg:10.2f}  REU")
    print(f"  buried SASA        {sasa:10.1f}  A^2")
    print(f"  dG per 1000 A^2    {1000*dg/sasa:10.2f}")
    print(f"  interface residues {ia.get_num_interface_residues():10d}")
    print(f"  packstat           {d.packstat:10.3f}")
    print(f"  H-bond energy      {d.total_hb_E:10.2f}")

analyse('ppk2_pair_CD.pdb', 'C-D  (923 A2, smallest)')
