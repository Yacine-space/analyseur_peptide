from services import analyses
from config import AA_DATABASE
from visualisation import plot

class Peptide:
    def __init__(self, sequence, n_terminal, c_terminal):
        self.sequence = sequence.strip().upper()
        self.n_terminal = n_terminal
        self.c_terminal = c_terminal
        self._valider()

    def _valider(self):
        invalides = {aa for aa in self.sequence if aa not in AA_DATABASE}
        if invalides:
            raise ValueError(f"Acides aminés invalides : {','.join(sorted(invalides))}")

    @property
    def masse_molaire(self):
        return analyses.masse_mol(self)
    
    @property
    def charge_nette_ph_7(self):
        return analyses.charge_nette_ph_7(self)

    @property
    def letter_code(self):
        return analyses.letter_code(self)

    @property
    def coefficient_extinction(self):
        return analyses.coefficient_extinction(self)

    @property
    def solubilite(self):
        return analyses.solubilite(self)

    @property
    def hydrophilie_moyenne(self):
        return analyses.hydrophilie_moyenne(self)
    
    @property
    def longueur(self):
        return len(self.sequence)

    @property
    def phi(self):
        return analyses.calc_phi(self)

    def get_courbe_charge(self):
        return plot.graphique_charge_ph(self)
    
    def get_courbe_hydrophilicite(self):
        return plot.graphique_hydrophilicite(self)
