import joblib
import numpy as np
from models.peptide import Peptide
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class AMPPredictor:

    def __init__(self):
        self.model = joblib.load(os.path.join(BASE_DIR, "data", "model.pkl"))
        self.scaler = joblib.load(os.path.join(BASE_DIR, "data", "scaler.pkl"))

    def _build_feature_vector(self, sequence):
        AA_ORDER = sorted("ACDEFGHIKLMNPQRSTVWY")
        p= Peptide(sequence, "no_modif", "no_modif")
        #physico-chimique
        physico = [
            p.masse_molaire,
            p.charge_nette_ph_7,
            p.hydrophilie_moyenne,
            p.phi,
            p.longueur,
            p.coefficient_extinction
        ]
        nterm_seq= sequence[:15].ljust(15, "-")
        binary_nterm=[]
        for a in nterm_seq:
            for aa in AA_ORDER:
                if a == aa:
                    binary_nterm.append(1)
                else:
                    binary_nterm.append(0)
            
        cterm_seq = sequence[-15:].ljust(15, "-")
        binary_cterm=[]
        for a in cterm_seq:
            for aa in AA_ORDER:
                if a == aa:
                    binary_cterm.append(1)
                else:
                    binary_cterm.append(0)

        return physico + binary_nterm + binary_cterm
    
    def predire(self, peptide: Peptide) -> dict:
        features = self._build_feature_vector(peptide.sequence)
        features_scaled = self.scaler.transform([features])

        label = self.model.predict(features_scaled)[0]
        probabilite = self.model.predict_proba(features_scaled)[0][1]

        return {
            "prediction" : "AMP" if label == 1 else "Non-AMP",
            "probabilite" : round(float(probabilite) * 100, 1) 
        }
        
