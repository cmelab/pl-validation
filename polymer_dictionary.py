# Dictionary of polymers and their key parameters

bond_indices = {
    "PCPDTPT_nC16": {"Bond indices": [121,122], "Orientations":[[1,0,0],[-1,0,0]], "Density": 0.05},
    "PCPDTFBT_C1_BO": {"Bond indices": [196,197], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05},
    "PCPDTFBT_C5_BO": {"Bond indices": [244,245], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05},
    "PCPDT_PT_eneHD": {"Bond indices": [122,123], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.005},
    "PCPDTPT_eneODD": {"Bond indices": [12,147], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05},
    "PCPDTFBT_C3_BO": {"Bond indices": [220,221], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05},
    "PCPDTFBT_C11_BO": {"Bond indices": [316,317], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.005},
    "PIDTBT_nC16": {"Bond indices": [229,230], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05},
    "PIDTFBT_C11_BO": {"Bond indices": [602,603], "Orientations":[[1,0,0],[-1,0,0]], "Density": 0.005},
    "PIDTCPDT_C11BO": {"Bond indices": [437,438], "Orientations":[[-1,0,0],[1,0,0]],"Density": 0.05},
    "PCPDTPT_HD": {"Bond indices": [121,122], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05},
    "PCPDTPT_ODD": {"Bond indices": [145,146], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05},
    "PCPDTFBT_C4_BO": {"Bond indices": [232,233], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.01}
}

molecules_not_tested_in_paper = {
    "PDPPPyT_ODD": {"Bond indicies": [164,165], "Orientations":[[1,1,1],[-1,-1,-1]]},
    "fixed_PCPDTPT_nC16": {"Bond indices": [52,53], "Orientations":[[-1,0,0],[1,0,0]], "Density": "0.05"} # Nitrogen in the correct position
}

Fixed_PCPDTPT_nC16_SMILES_string = "C1SC2C3SC(C4C5=NSN=C5C=NC=4)=CC=3C(CCCCCCCCCCCCCCCC)(CCCCCCCCCCCCCCCC)C=2C=1" # Fairly certain the carbon rings should be aromatic
aromatic_version = "c1Sc2c3Sc(c4c5=NSN=c5c=Nc=4)=cc=3c(CCCCCCCCCCCCCCCC)(CCCCCCCCCCCCCCCC)c=2c=1" # Running into RDKit error when parameterizing