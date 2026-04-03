from config import AA_DATABASE
import numpy as np
def ph(peptide):
    ph = -1
    d={}
    while ph<=13:
        ph+=1
        for aa in peptide.sequence:
            charge=infos_peptid(peptide, ph)
            d[ph]=charge
    return d

def infos_peptid(peptide, ph):
    charge_N_term=0
    charge_C_term=0
    fra_base=0
    fra_acide=0
    for aa in peptide.sequence:
        if aa in AA_DATABASE:
            if AA_DATABASE[aa]["type"]=="acide":
                pka = AA_DATABASE[aa]["pKa"]
                charge_acid= 1/(1+10**(pka-ph))
                fra_acide+=charge_acid
            elif AA_DATABASE[aa]["type"]=="base":
                pka= AA_DATABASE[aa]['pKa']
                charge_base=1/(1+10**(ph-pka))
                fra_base+=charge_base

        charge_nette=fra_base-fra_acide
        charge_nette+= -1/(1+10**(2-ph)) #extrimité c-ter
        charge_nette+= 1/(1+10**(ph-9))  #extrimité N--ter
    charge_N_term, charge_C_term,_,_ = modification(peptide)   #ajout charges des modifications N et C term
    charge_nette=charge_nette+charge_N_term+charge_C_term 
    return charge_nette

def masse_mol(peptide):
    _,_,poids_mol_N,poids_mol_C=modification(peptide)
    masse_mol_tot=0
    for aa in peptide.sequence:
        if aa in AA_DATABASE:
            masse_mol_tot+=AA_DATABASE[aa]["masse_moleculaire"]
    masse_mol_pep=masse_mol_tot-(len(peptide.sequence)-1)*18.015 #soustraire les liaisons hydro
    masse_mol_pep=masse_mol_pep+poids_mol_N+poids_mol_C #ajout masse moléculaire des modifications N et C term
    masse_mol_pep=round(masse_mol_pep, 2)
    return masse_mol_pep

def calc_phi(peptide):
    ph=0
    e={}
    while ph<14:
        ph+=0.01
        ph= round(ph, 2)
        phi=infos_peptid(peptide, ph)
        e[ph]=phi
        isoelectric_point=0
        peite_valeur_charge=float('inf')
        for key, value in e.items():
            if value>0:
                if value<peite_valeur_charge:
                    isoelectric_point=key
    return  isoelectric_point
            
def hydrophilie_moyenne(peptide):
    hydrophilicite_tot = []
    for aa in peptide.sequence:
        if aa in AA_DATABASE:
            hydrophilicite_tot.append(AA_DATABASE[aa]['hydrophilicte'])

        tot= sum(hydrophilicite_tot)
        hydro_moy= tot/len(peptide.sequence) 
        hydro_moy= round(hydro_moy,2)
    return hydro_moy

def charge_nette_ph_7(peptide):
    charge_nette= ph(peptide)
    charge_nette_ph_7= charge_nette[7]
    charge_nette_ph_7= round(charge_nette_ph_7, 2)
    return charge_nette_ph_7

    
def letter_code(peptide):
    letter_code=[]
    for aa in peptide.sequence:
        if aa in AA_DATABASE:
            letter_code.append(AA_DATABASE[aa]["code"])
    letter_code_str="-".join(letter_code)
    return letter_code_str

def charge_vs_ph(peptide):
    valeurs_ph= np.arange(0, 14.1, 0.1)
    charges = [infos_peptid(peptide, ph) for ph in valeurs_ph]
    return valeurs_ph, charges

def aa_vs_hydrophilicte (peptide):
    aa_list= []
    hydro_list= []
    couleur = []
    for aa in peptide.sequence:
        aa_list.append(aa)
        hydro_list.append(AA_DATABASE[aa]['hydrophilicte'])
    
    for v in hydro_list:
        if v < 0:
            couleur.append('red')
        elif v >0:
            couleur.append('green')
        else:
            couleur.append('grey')
    
    return aa_list, hydro_list, couleur

def modification(peptide):
    modificatons_n = {
        "no_modif"     : (0, 0),
        "acetylation"  : (-1, 43),
        "biotine"      : (0, 28),
        "formylation"  : (0, 28),
        "dansyl"       : (0, 249),
        "pyroglutamate": (-1, -18),
    }
    modifications_c = {
        "no_modif"     : (0, 0),
        "amidation"    : (1, -17),
        "cyclisation"  : (0, -18),
        "esterification": (0, 14),
        "biotine"      : (0, 340),
    }
    charge_N_term, poids_mol_N = modificatons_n.get(peptide.n_terminal, (0, 0))
    charge_C_term, poids_mol_C = modifications_c.get(peptide.c_terminal, (0, 0))
    
    return charge_N_term, charge_C_term, poids_mol_N, poids_mol_C

def coefficient_extinction(peptide):
    coefficient_extinction=0
    for aa in peptide.sequence:
        if aa == "W":
            coefficient_extinction+=5690
        elif aa == "Y":
            coefficient_extinction+=1280
    return coefficient_extinction

def score_solubilite(peptide):
    aa_hydrophobes = 0
    aa_polaires = 0
    hydrophobes = set(['A', 'V', 'L', 'I', 'M', 'F', 'W'])
    polaires = set(['S', 'T', 'N', 'Q', 'K', 'R', 'D', 'E'])
    for aa in peptide.sequence:
            if aa in hydrophobes:
                aa_hydrophobes+=1
            elif aa in polaires:
                aa_polaires+=1
    longueur_pep = len(peptide.sequence)
    pourc_hydrophobes = (aa_hydrophobes/longueur_pep)*100
    pourc_polaires = (aa_polaires/longueur_pep)*100
    charge_nette= charge_nette_ph_7(peptide)
    score_solu=(charge_nette*5)-pourc_hydrophobes+pourc_polaires
    return score_solu

def solubilite(peptide):
    score_solu=score_solubilite(peptide)
    if score_solu>100:
        return "Excellente solubilité"
    elif 50<score_solu<100:
        return "Bonne solubilité"
    elif 0<score_solu<50:
        return "Solubilité modérée"
    elif score_solu<0:
        return "Faible solubilité (risque d'agrégation)"