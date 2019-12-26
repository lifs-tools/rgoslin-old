from pygoslin.parser.BaseParserEventHandler import BaseParserEventHandler
from pygoslin.domain.LipidAdduct import LipidAdduct
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.Adduct import Adduct
from pygoslin.domain.MolecularFattyAcid import MolecularFattyAcid
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo
from pygoslin.domain.LipidSpecies import LipidSpecies
from pygoslin.domain.LipidMolecularSubspecies import LipidMolecularSubspecies
from pygoslin.domain.LipidStructuralSubspecies import LipidStructuralSubspecies
from pygoslin.domain.StructuralFattyAcid import StructuralFattyAcid

class GoslinParserEventHandler(BaseParserEventHandler):
    
    def __init__(self):
        super().__init__()
        self.reset_lipid(None)
        
        self.registered_events["lipid_pre_event"] = self.reset_lipid
        self.registered_events["lipid_post_event"] = self.build_lipid
        
        self.registered_events["hg_cl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_mlcl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_pl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_lpl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_lpl_o_pre_event"] = self.set_head_group_name
        self.registered_events["hg_pl_o_pre_event"] = self.set_head_group_name
        self.registered_events["hg_lsl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_dsl_pre_event"] = self.set_head_group_name
        self.registered_events["ch_pre_event"] = self.set_head_group_name
        self.registered_events["hg_che_pre_event"] = self.set_head_group_name
        self.registered_events["mediator_pre_event"] = self.set_head_group_name
        self.registered_events["hg_mgl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_dgl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_sgl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_tgl_pre_event"] = self.set_head_group_name
        
        self.registered_events["gl_species_pre_event"] = self.set_species_level
        self.registered_events["pl_species_pre_event"] = self.set_species_level
        self.registered_events["chc_pre_event"] = self.set_species_level
        self.registered_events["sl_species_pre_event"] = self.set_species_level
        self.registered_events["fa2_unsorted_pre_event"] = self.set_molecular_subspecies_level
        self.registered_events["fa3_unsorted_pre_event"] = self.set_molecular_subspecies_level
        self.registered_events["fa4_unsorted_pre_event"] = self.set_molecular_subspecies_level
        
        self.registered_events["lcb_pre_event"] = self.new_lcb
        self.registered_events["lcb_post_event"] = self.clean_lcb
        self.registered_events["fa_pre_event"] = self.new_fa
        self.registered_events["fa_post_event"] = self.append_fa
        
        self.registered_events["ether_pre_event"] = self.add_ether
        self.registered_events["old_hydroxyl_pre_event"] = self.add_old_hydroxyl
        self.registered_events["db_count_pre_event"] = self.add_double_bonds
        self.registered_events["carbon_pre_event"] = self.add_carbon
        self.registered_events["hydroxyl_pre_event"] = self.add_hydroxyl
        
        
        self.registered_events["adduct_info_pre_event"] = self.new_adduct
        self.registered_events["adduct_pre_event"] = self.add_adduct
        self.registered_events["charge_pre_event"] = self.add_charge
        self.registered_events["charge_sign_pre_event"] = self.add_charge_sign
    

    def reset_lipid(self, node):
        self.level = LipidLevel.STRUCTURAL_SUBSPECIES
        self.lipid = None
        self.head_group = ""
        self.lcb = None
        self.fa_list = []
        self.current_fa = None
        self.adduct = None
        

    def set_head_group_name(self, node):
        self.head_group = node.get_text()
        
    
    def set_species_level(self, node):
        self.level = LipidLevel.SPECIES
        
        
    def set_molecular_subspecies_level(self, node):
        self.level = LipidLevel.MOLECULAR_SUBSPECIES
        
        
    def new_fa(self, node):
        if self.level == LipidLevel.SPECIES:
            self.current_fa = LipidSpeciesInfo()
            
        elif self.level == LipidLevel.MOLECULAR_SUBSPECIES:
            self.current_fa = MolecularFattyAcid("FA%i" % (len(self.fa_list) + 1), 2, 0, 0, LipidFaBondType.ESTER, False, -1)
            
        elif self.level == LipidLevel.STRUCTURAL_SUBSPECIES:
            self.current_fa = StructuralFattyAcid("FA%i" % (len(self.fa_list) + 1), 2, 0, 0, LipidFaBondType.ESTER, False, 0)
        
        
    def new_lcb(self, node):
        if self.level == LipidLevel.SPECIES:
            self.lcb = LipidSpeciesInfo()
            self.lcb.lipid_FA_bond_type = LipidFaBondType.ESTER
            
        elif self.level == LipidLevel.STRUCTURAL_SUBSPECIES:
            self.lcb = StructuralFattyAcid("LCB", 2, 0, 1, LipidFaBondType.ESTER, True, 1)
            
        self.current_fa = self.lcb
            
            
    def clean_lcb(self, node):
        self.current_fa = None
        
        
            
    def append_fa(self, node):
        if self.level == LipidLevel.STRUCTURAL_SUBSPECIES:
            self.current_fa.position = len(self.fa_list) + 1
            
        self.fa_list.append(self.current_fa)
        self.current_fa = None
        
        
    def build_lipid(self, node):
        if self.lcb != None:
            for fa in self.fa_list: fa.position += 1
            self.fa_list = [self.lcb] + self.fa_list
        
        lipid = None
        
        
        
        if self.level == LipidLevel.SPECIES:
            if len(self.fa_list) > 0:
                lipid_species_info = LipidSpeciesInfo(self.fa_list[0])
                lipid_species_info.level = LipidLevel.SPECIES
                lipid = LipidSpecies(self.head_group, lipid_species_info = lipid_species_info)
            else:
                lipid = LipidSpecies(self.head_group)
            
        elif self.level == LipidLevel.MOLECULAR_SUBSPECIES:
            lipid = LipidMolecularSubspecies(self.head_group, self.fa_list)
            
        elif self.level == LipidLevel.STRUCTURAL_SUBSPECIES:
            lipid = LipidStructuralSubspecies(self.head_group, self.fa_list)
    
        self.lipid = LipidAdduct()
        self.lipid.lipid = lipid
        self.lipid.adduct = self.adduct
        
        
    def add_ether(self, node):
        ether = node.get_text()
        if ether == "a": self.current_fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMANYL
        elif ether == "p": self.current_fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMENYL
        
        
    def add_old_hydroxyl(self, node):
        old_hydroxyl = node.get_text()
        if old_hydroxyl == "d": self.current_fa.num_hydroxyl = 2
        if old_hydroxyl == "t": self.current_fa.num_hydroxyl = 3
        
        
    def add_double_bonds(self, node):
        self.current_fa.num_double_bonds = int(node.get_text())
        
        
    def add_carbon(self, node):
        self.current_fa.num_carbon = int(node.get_text())
        
        
    def add_hydroxyl(self, node):
        self.current_fa.num_hydroxyl = int(node.get_text())
        
        
    def new_adduct(self, node):
        self.adduct = Adduct("", "", 0, 0)
        
        
    def add_adduct(self, node):
        self.adduct.adduct_string = node.get_text()
        
        
    def add_charge(self, node):
        self.adduct.charge = int (node.get_text())
        
        
    def add_charge_sign(self, node):
        sign = node.get_text()
        if sign == "+": self.adduct.set_charge_sign(1)
        if sign == "-": self.adduct.set_charge_sign(-1)
        
        
        