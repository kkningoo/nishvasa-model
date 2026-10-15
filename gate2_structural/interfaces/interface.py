from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.core.scoring import get_score_function

init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false")
sf = get_score_function()

def analyse(pdb, label):
    pose = pose_from_pdb(pdb)
    # jump 1 separates the two chains
    ia = InterfaceAnalyzerMover(1, False, sf, True, False, True, False)
    ia.apply(pose)
    d = ia.get_all_data()
    dg   = ia.get_interface_dG()
    sasa = ia.get_interface_delta_sasa()
    print(f"\n=== {label} ===")
    print(f"  dG_separated       {dg:10.2f}  REU")
    print(f"  buried SASA        {sasa:10.1f}  A^2")
    print(f"  dG per 1000 A^2    {1000*dg/sasa:10.2f}")
    print(f"  interface residues {ia.get_num_interface_residues():10d}")
    print(f"  packstat           {d.packstat:10.3f}")
    print(f"  H-bond energy      {d.total_hb_E:10.2f}")

analyse('gs_pair_AB.pdb', 'intra-ring A-B')
analyse('gs_pair_AH.pdb', 'CROSS-RING A-H')
