# Dictionary of polymers and their key parameters
import ele

polymer_dictionary = {
    "PCPDTFBT_C11_BO": {"Bond indices": [316,317], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.005, "Temperatures":[0,0,0], "Sulfur index":28},
    "PCPDTFBT_C1_BO": {"Bond indices": [196,197], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05, "Temperatures":[0,0,0], "Sulfur index":28},
    "PCPDTFBT_C3_BO": {"Bond indices": [220,221], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05, "Temperatures":[0,0,0], "Sulfur index":28},
    "PCPDTFBT_C4_BO": {"Bond indices": [232,233], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.01, "Temperatures":[0,0,0], "Sulfur index":28},
    "PCPDTFBT_C5_BO": {"Bond indices": [244,245], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05, "Temperatures":[0,0,0], "Sulfur index":28},    
    "PCPDTPT_HD": {"Bond indices": [121,122], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05, "Temperatures":[0,0,0], "Sulfur index":10},
    "PCPDTPT_ODD": {"Bond indices": [145,146], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05, "Temperatures":[0,0,0], "Sulfur index":10},
    "PCPDTPT_eneODD": {"Bond indices": [12,147], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05, "Temperatures":[0,0,0], "Sulfur index":10},
    "PCPDTPT_nC16": {"Bond indices": [121,122], "Orientations":[[1,0,0],[-1,0,0]], "Density": 0.05, "Temperatures":[0,0,0], "Sulfur index":10},
    "PCPDT_PT_eneHD": {"Bond indices": [122,123], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.005, "Temperatures":[0,0,0], "Sulfur index":16},
    "PIDTBT_nC16": {"Bond indices": [229,230], "Orientations":[[-1,0,0],[1,0,0]], "Density": 0.05, "Temperatures":[0,0,0], "Sulfur index":15},
    "PIDTCPDT_C11BO": {"Bond indices": [437,438], "Orientations":[[-1,0,0],[1,0,0]],"Density": 0.05, "Temperatures":[0,0,0], "Sulfur index":15},
    "PIDTFBT_C11_BO": {"Bond indices": [602,603], "Orientations":[[1,0,0],[-1,0,0]], "Density": 0.005, "Temperatures":[0,0,0], "Sulfur index":15}
}

molecules_not_tested_in_paper = {
    "PDPPPyT_ODD": {"Bond indicies": [164,165], "Orientations":[[1,1,1],[-1,-1,-1]]},
    "alternate_PCPDTPT_nC16": {"Bond indices": [52,53], "Orientations":[[-1,0,0],[1,0,0]], "Density": "0.05"} # Nitrogen in the correct position
}

alternate_PCPDTPT_nC16_SMILES_string = "C1SC2C3SC(C4C5=NSN=C5C=NC=4)=CC=3C(CCCCCCCCCCCCCCCC)(CCCCCCCCCCCCCCCC)C=2C=1" # P and PT carbon rings should be aromatic
aromatic_version = "c1Sc2c3Sc(c4c5=NSN=c5c=Nc=4)=cc=3c(CCCCCCCCCCCCCCCC)(CCCCCCCCCCCCCCCC)c=2c=1" # Not correct aromaticity

polymer_fragments = {
    "P":"c1ccncc1",
    "CPDT":"C1C3=C(SC=C3)C2=C1C=CS2",
    "PT":"n1ccc2nsnc2c1",
    "IDT":"C1C2=C(C3=CC4=C(C=C31)C5=C(C4)C=CS5)SC=C2",
    "BT":"c1ccc2nsnc2c1",
    "DPP":"N1C=C(C(=O)Nc2)c2C1(=O)",
    "FBT":"c1(F)ccc2nsnc2c1",
    "T":"c1ccsc1",
    "CPDT_eneHD":"C1(=C(CC(CCCCCC)CCCCCCCC)CC(CCCCCC)CCCCCCCC)C3=C(SC=C3)C2=C1C=CS2",
    "TPD":"C1(=O)C2=CSC=C2C(=O)N1(CCCCCCCC)",
    "BDT":"s1c2c(OCC(CC)CCCC)c3ccsc3c(OCC(CC)CCCC)c2cc1"
}

element_dict = {("C0"): ele.element_from_symbol("C"),
                 ("C1"): ele.element_from_symbol("C"),
                 ("C3"): ele.element_from_symbol("C"),
                 ("C5"): ele.element_from_symbol("C"),
                 ("C6"): ele.element_from_symbol("C"),
                 ("N2"): ele.element_from_symbol("N"),
                 ("N3"): ele.element_from_symbol("N"),
                 ("N5"): ele.element_from_symbol("N"),
                 ("F4"): ele.element_from_symbol("F"),
                 ("S1"):ele.element_from_symbol("S"),
                 ("S2"):ele.element_from_symbol("S")}