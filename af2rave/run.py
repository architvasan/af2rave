import af2rave.alphafold as alphafoldrv
import af2rave.feature as featurerv
import af2rave.amino as aminorv
import af2rave.simulation as simulationrv
import af2rave.spib as spibrv
from dataclasses import dataclass
import glob 
"""
## Class to run AF2RAVE 
### Parameters:
- sequence: str, protein target sequence 
- jobname: str 
- output_af2dir: str, output for af2
-homooligomer: str 
- ref_pdb: str, reference for rmsd filtering 
- rmsd_sel: str, selection for rmsd 
- rmsd_cutoff: how much rmsd to filter by 
"""
@dataclass
class RunAF2RAVE:
    sequence: str 
    jobname: str 
    output_af2dir: str 
    homooligomer: str = "1" 
    ref_pdb: object=None
    rmsd_sel: str = "name CA"
    rmsd_cutoff: float = 10.0
    cluster_sel: str | tuple[str, str] = "name CA"
    """
    Run colabfold with reduced msa to get conformations
    """
    def run_colabfold(self):
        cf_fold_cl = alphafoldrv.colabfold.ColabFold(
                                    sequence=self.sequence,
                                    name = self.jobname,
                                    output_dir = self.output_af2dir)
        cf_fold_cl.predict()

    """
    apply rmsd based filter to filter structures that are irrelevant
    """
    def run_rmsd_filter(self):
        self.pdb_list = glob.glob(f'{self.output_af2dir}/*pdb') 
        self.feature_selection_cl = featurerv.analysis.FeatureSelection(
                                   self.pdb_list,
                                   self.ref_pdb                 
                                )
        mask = self.feature_selection_cl.rmsd_filter(
                                        selection=self.rmsd_sel,
                                        rmsd_cutoff=self.rmsd_cutoff)
        self.feature_selection_cl.apply_filter(mask)
    
    """
    Cluster based on a selection (use a priori information)
    """
    def run_clustering(self):
        if hasattr(self, 'feature_selection_cl'):
            pass 
        else:
            self.pdb_list = glob.glob(f'{self.output_af2dir}/*pdb') 
            self.feature_selection_cl = featurerv.analysis.FeatureSelection(
                                   self.pdb_list,
                                   self.ref_pdb                 
                                ) 
        names_sorted, coeffvars = self.feature_selection_cl.rank_feature(
                                        selection = self.cluster_sel,
                                               )
        center_ids = self.feature_selection_cl.regular_space_clustering(
                                                names_sorted,
                                                min_dist = 0.5
                                        )
        self.pdb_centers = [self.feature_selection_cl.pdb_name[cit]\
                                for cit in center_ids]

    def run_simulation(self,
                       pdb_cent_list: list[str]| None = None):
        if pdb_cent_list ==None:
            pdb_cent_list = self.pdb_centers
        for pdb_it in pdb_cent_list:
            simbox_cl = simulationrv.SimulationBox(pdb_it)
            simbox_cl.create_box()
            simbox_cl.save_pdb(f'{pdb_it}.updated.pdb')
            sim_cl = simulationrv.UnbiasedSimulation(f'{pdb_it}.updated.pdb')
            sim_cl.run(5000)


if __name__ == "__main__":
    sequence = 'MQRGKVKWFNNEKGYGFIEVEGGSDVFVHFTAIQGEGFKTLEEGQEVSFEIVQGNRGPQAANVVKE' #@param {type:"string"}
    jobname = 'CSP'
    output_dir = 'output_test'
    runner = RunAF2RAVE(
                sequence,
                jobname,
                output_dir,
    )
    runner.run_rmsd_filter()
    runner.run_clustering()
    runner.run_simulation()
    
    print(len(runner.pdb_centers))