import os
import pickle
from scipy.constants import e
import molecular_function as mf

asset_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'asset')
model_dir = os.path.join(asset_dir, 'model')


def ml_dock(data):
    # Load model
    with open(os.path.join(model_dir, 'Linear_Regression_Model.pkl'), "rb") as f:
        model = pickle.load(f)

    # Make predictions
    predictions = model.predict([data])

    return predictions[0]

class docking:
    def __init__(self, ligand="", target=""):
        self.ligand = ligand
        self.target = target

    def run(self):
        ligand = self.ligand
        target = self.target
        aa1 = mf.molecular_function.calculate_amino_acid_center_of_mass(ligand)
        smiles1 = mf.molecular_function.seq_to_smiles(ligand)
        aa2 = mf.molecular_function.calculate_amino_acid_center_of_mass(target)
        smiles2 = mf.molecular_function.seq_to_smiles(target)
        dist1 = mf.molecular_function.calculate_distance_between_amino_acids(aa1, aa2)
        attract = mf.ms.attractive_energy(dist1)
        repulsive = mf.ms.repulsive_energy(dist1)
        vdw = mf.ms.lj_force(dist1)
        ce = mf.ms.coulomb_energy(e, e, dist1)
        ff = vdw+ce
        molwt1 = mf.molecular_function.calculate_molecular_weight(smiles1)
        molwt2 = mf.molecular_function.calculate_molecular_weight(smiles2)
        dist2 = mf.molecular_function.calculate_distance_between_amino_acids(aa1, aa2)
        feature = [aa1, molwt1, aa2, molwt2, dist2]
        mldock_res = ml_dock(feature)

        res = {
            'force_field' : ff,
            'ml_dock' : mldock_res
        }

        return res