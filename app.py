import matplotlib.pyplot as plt
import numpy as np
from flask import Flask, request, render_template
import plotly.graph_objs as go
import plotly.io as pio
import os

proprietes_aa = {
    "A": {"nom": "Alanine","code": "ALA","masse_moleculaire": 89.09, "hydrophilicte": -0.5, "type": "neutre"},
    "R": {"nom": "Arginine","code": "ARG","masse_moleculaire": 174.20, "hydrophilicte": 3.0, "type": "base", "pKa": 12.5},
    "N": {"nom": "Asparagine","code": "ASN","masse_moleculaire": 132.12, "hydrophilicte": 0.2, "type": "neutre"},
    "D": {"nom": "Aspartate","code": "ASP","masse_moleculaire": 133.10, "hydrophilicte": 3.0, "type": "acide", "pKa": 3.9},
    "C": {"nom": "Cysteine","code": "CYS","masse_moleculaire": 121.15, "hydrophilicte": -1.0, "type": "acide", "pKa": 8.3},
    "Q": {"nom": "Glutamine","code": "GLN","masse_moleculaire": 146.15, "hydrophilicte": 0.2, "type": "neutre"},
    "E": {"nom": "Glutamate","code": "GLU","masse_moleculaire": 147.13, "hydrophilicte": 3.0, "type": "acide", "pKa": 4.1},
    "G": {"nom": "Glycine","code": "GLY","masse_moleculaire": 75.07, "hydrophilicte": 0.0, "type": "neutre"},
    "H": {"nom": "Histidine","code": "HIS","masse_moleculaire": 155.16, "hydrophilicte": -0.5, "type": "base", "pKa": 6.0},
    "I": {"nom": "Isoleucine","code": "ILE","masse_moleculaire": 131.17, "hydrophilicte": -1.8, "type": "neutre"},
    "L": {"nom": "Leucine","code": "LEU","masse_moleculaire": 131.17, "hydrophilicte": -1.8, "type": "neutre"},
    "K": {"nom": "Lysine","code": "LYS","masse_moleculaire": 146.19, "hydrophilicte": 3.0, "type": "base", "pKa": 10.5},
    "M": {"nom": "Methionine","code": "MET","masse_moleculaire": 149.21, "hydrophilicte": -1.3, "type": "neutre"},
    "F": {"nom": "Phenylalanine","code": "PHE","masse_moleculaire": 165.19, "hydrophilicte": -2.5, "type": "neutre"},
    "P": {"nom": "Proline","code": "PRO","masse_moleculaire": 115.13, "hydrophilicte": 0.0, "type": "neutre"},
    "S": {"nom": "Serine","code": "SER","masse_moleculaire": 105.09, "hydrophilicte": 0.3, "type": "neutre"},
    "T": {"nom": "Threonine","code": "THR","masse_moleculaire": 119.12, "hydrophilicte": -0.4, "type": "neutre"},
    "W": {"nom": "Tryptophane","code": "TRP","masse_moleculaire": 204.23, "hydrophilicte": -3.4, "type": "neutre"},
    "Y": {"nom": "Tyrosine","code": "TYR","masse_moleculaire": 181.19, "hydrophilicte": -2.3, "type": "acide", "pKa": 10.1},
    "V": {"nom": "Valine","code": "VAL","masse_moleculaire": 117.15, "hydrophilicte": -1.5, "type": "neutre"},
}

extremites = {
    "N_terminale": {"pKa": 9.0, "type": "base"},
    "C_terminale": {"pKa": 2.0, "type": "acide"}
}

def ph(seq_peptide, proprietes_aa, N_terminal, C_terminal):
    ph = -1
    d={}
    
    while ph<=13:
        ph+=1
        for aa in seq_peptide:
            charge=infos_peptid(seq_peptide, proprietes_aa, ph, N_terminal, C_terminal)
            d[ph]=charge

    return d

def infos_peptid(seq_peptide, proprietes_aa,ph, N_terminal, C_terminal):
    charge_N_term=0
    charge_C_term=0
    fra_base=0
    fra_acide=0
    for aa in seq_peptide:
        if aa in proprietes_aa:
            if proprietes_aa[aa]["type"]=="acide":
                pka = proprietes_aa[aa]["pKa"]
                charge_acid= 1/(1+10**(pka-ph))
                fra_acide+=charge_acid
            elif proprietes_aa[aa]["type"]=="base":
                pka= proprietes_aa[aa]['pKa']
                charge_base=1/(1+10**(ph-pka))
                fra_base+=charge_base

        charge_nette=fra_base-fra_acide
        charge_nette+= -1/(1+10**(2-ph)) #extrimité c-ter
        charge_nette+= 1/(1+10**(ph-9))  #extrimité N--ter
    charge_N_term, charge_C_term,_,_ = modification(N_terminal, C_terminal)   #ajout charges des modifications N et C term
    charge_nette=charge_nette+charge_N_term+charge_C_term 
    return charge_nette

def masse_mol(seq_peptide, proprietes_aa, N_terminal, C_terminal):
    _,_,poids_mol_N,poids_mol_C=modification(N_terminal,C_terminal)
    masse_mol_tot=0
    for aa in seq_peptide:
        if aa in proprietes_aa:
            masse_mol_tot+=proprietes_aa[aa]["masse_moleculaire"]

    masse_mol_pep=masse_mol_tot-(len(seq_peptide)-1)*18.015 #soustraire les liaisons hydro
    masse_mol_pep=masse_mol_pep+poids_mol_N+poids_mol_C #ajout masse moléculaire des modifications N et C term
    masse_mol_pep=round(masse_mol_pep, 2)
    return masse_mol_pep

def calc_phi(seq_peptide, proprietes_aa, N_terminal, C_terminal):
    ph=0
    e={}
    while ph<14:
        ph+=0.01
        ph= round(ph, 2)
        phi=infos_peptid(seq_peptide, proprietes_aa, ph, N_terminal, C_terminal)
        e[ph]=phi
        isoelectric_point=0
        peite_valeur_charge=float('inf')
        for key, value in e.items():
            if value>0:
                if value<peite_valeur_charge:
                    isoelectric_point=key
    return  isoelectric_point
            
def hydrophilie_moyenne(seq_peptide, proprietes_aa):
    hydrophilicite_tot = []
    for aa in seq_peptide:
        if aa in proprietes_aa:
            hydrophilicite_tot.append(proprietes_aa[aa]['hydrophilicte'])

        tot= sum(hydrophilicite_tot)
        hydro_moy= tot/len(seq_peptide) 
        hydro_moy= round(hydro_moy,2)
    return hydro_moy

def charge_nette_ph_7(seq_peptide, proprietes_aa, N_terminal, C_terminal):
    charge_nette= ph(seq_peptide, proprietes_aa, N_terminal, C_terminal)
    charge_nette_ph_7= charge_nette[7]
    charge_nette_ph_7= round(charge_nette_ph_7, 2)
    return charge_nette_ph_7

def length_pep(seq_peptide):
    length_pep=len(seq_peptide)
    return length_pep
    
def letter_code(seq_peptide, proprietes_aa):
    letter_code=[]
    for aa in seq_peptide:
        if aa in proprietes_aa:
            letter_code.append(proprietes_aa[aa]["code"])
    letter_code_str="-".join(letter_code)
    return letter_code_str

def charge_vs_ph(seq_peptide, proprietes_aa, N_terminal, C_terminal):
    valeurs_ph= np.arange(0, 14.1, 0.1)
    charges = [infos_peptid(seq_peptide, proprietes_aa, ph, N_terminal, C_terminal) for ph in valeurs_ph]
    return valeurs_ph, charges

def aa_vs_hydrophilicte (seq_peptide, proprietes_aa):
    aa_list= []
    hydro_list= []
    couleur = []
    for aa in seq_peptide:
        aa_list.append(aa)
        hydro_list.append(proprietes_aa[aa]['hydrophilicte'])
    
    for v in hydro_list:
        if v < 0:
            couleur.append('red')
        elif v >0:
            couleur.append('green')
        else:
            couleur.append('grey')
    
    return aa_list, hydro_list, couleur

def modification(N_terminal, C_terminal):
    if N_terminal == "no_modif":
        charge_N_term=0
        poids_mol_N=0
    elif N_terminal == "acetylation":
        charge_N_term= -1
        poids_mol_N= 43
    elif N_terminal == "biotine":
        charge_N_term= 0
        poids_mol_N=28
    elif N_terminal=="formylation":
        charge_N_term= 0
        poids_mol_N=28
    elif N_terminal=="dansyl":
        charge_N_term= 0
        poids_mol_N=249
    elif N_terminal=="pyroglutamate":
        charge_N_term=-1
        poids_mol_N=-18
    if C_terminal=="no_modif":
        charge_C_term=0
        poids_mol_C=0
    elif C_terminal== "amidation":
        charge_C_term= 1
        poids_mol_C= -17
    elif C_terminal == "cyclisation":
        charge_C_term= 0
        poids_mol_C= -18
    elif C_terminal=="esterification":
        charge_C_term= 0
        poids_mol_C= 14
    elif C_terminal=="biotine":
        charge_C_term= 0
        poids_mol_C= 340
    return charge_N_term, charge_C_term, poids_mol_N, poids_mol_C

def coefficient_extinction(seq_peptide):
    coefficient_extinction=0
    for aa in seq_peptide:
        if aa == "W":
            coefficient_extinction+=5690
        elif aa == "Y":
            coefficient_extinction+=1280
    return coefficient_extinction

def score_solubilite(seq_peptide, N_terminal, C_terminal):
    aa_hydrophobes = 0
    aa_polaires = 0
    hydrophobes = set(['A', 'V', 'L', 'I', 'M', 'F', 'W'])
    polaires = set(['S', 'T', 'N', 'Q', 'K', 'R', 'D', 'E'])
    for aa in seq_peptide:
            if aa in hydrophobes:
                aa_hydrophobes+=1
            elif aa in polaires:
                aa_polaires+=1
    longueur_pep = len(seq_peptide)
    pourc_hydrophobes = (aa_hydrophobes/longueur_pep)*100
    pourc_polaires = (aa_polaires/longueur_pep)*100
    charge_nette= charge_nette_ph_7(seq_peptide, proprietes_aa, N_terminal, C_terminal)
    score_solu=(charge_nette*5)-pourc_hydrophobes+pourc_polaires
    return score_solu

def solubilite(seq_peptide, N_terminal, C_terminal):
    score_solu=score_solubilite(seq_peptide, N_terminal, C_terminal)
    if score_solu>100:
        return "Excellente solubilité"
    elif 50<score_solu<100:
        return "Bonne solubilité"
    elif 0<score_solu<50:
        return "Solubilité modérée"
    elif score_solu<0:
        return "Faible solubilité (risque d'agrégation)"

def graphique_hydrophilicite(seq_peptide):
    aa, hydrophilicte, couleur = aa_vs_hydrophilicte(seq_peptide, proprietes_aa)
    positions = list(range(1, len(seq_peptide) + 1))
    fig_hydro = go.Figure()
    fig_hydro.add_trace(go.Bar(
        x=positions,
        y=hydrophilicte,
        marker_color=couleur,
        text=aa,
        hovertemplate="AA: %{text}<br>Value: %{y}<extra></extra>"
    ))
    fig_hydro.update_layout(
        yaxis=dict(zeroline=True),
        xaxis=dict(
            dtick=1,
            tickmode="array",
            tickvals=positions,
            ticktext=aa,
        ),
        margin=dict(l=8, r=8, t=8, b=8),
        template="simple_white",
    )
    graph_2_html = pio.to_html(fig_hydro, full_html=False)
    return graph_2_html

def graphique_charge_ph(seq_peptide, N_terminal, C_terminal):
    valeurs_ph, charges = charge_vs_ph(seq_peptide, proprietes_aa, N_terminal, C_terminal)
    fig_charge = go.Figure()
    fig_charge.add_trace(go.Scatter(
        x=valeurs_ph,
        y=charges,
        mode="lines",
        line=dict(color="darkblue", width=3),
    ))

    fig_charge.update_layout(
        xaxis=dict(
            title={
                'text': 'pH',
                'font': {"size": 14, "color": "black"},
            },
            range=[0, 14],
            dtick=1,
            tickmode="linear",
            position=0.5,
            showline=True,
            linecolor="black",
        ),
        yaxis=dict(
            title={
                "text": "Charge nette",
                "font": {"size": 14, "color": 'black'},
            },
            range=[-6, 6],
            tickmode="linear",
            position=0.5,
            showline=True,
            linecolor="black",
        ),
        margin=dict(l=8, r=8, t=8, b=8),
        width=400,
        height=400
    )
    graph_1_html = pio.to_html(fig_charge, full_html=False)
    return graph_1_html
app = Flask(__name__)

@app.route('/', methods=["GET", "POST"])
def home():
    seq_peptide=""
    stat= None
    N_terminal=""
    C_terminal=""
    error_message = None
    graph_1_html= ""
    graph_2_html = ""
    if request.method == "POST":
        seq_peptide= request.form.get("peptid", "").strip().upper()
        N_terminal= request.form.get("N_term", "")
        C_terminal=request.form.get("C_term", "")
        if seq_peptide:
            invalides = sorted({aa for aa in seq_peptide if aa not in proprietes_aa})

            if invalides:
                error_message = (
                    "Les lettres suivantes ne correspondent a aucun acide amine valide : "
                    + ", ".join(invalides)
                )
            else:
                stat = {
                    "letter_code_1": seq_peptide,
                    'length': length_pep(seq_peptide),
                    "letter_code_3": letter_code(seq_peptide, proprietes_aa),
                    "poid_mol": masse_mol(seq_peptide, proprietes_aa, N_terminal, C_terminal),
                    "phi": calc_phi(seq_peptide, proprietes_aa, N_terminal, C_terminal),
                    "charge_nette_ph_7": charge_nette_ph_7(seq_peptide, proprietes_aa, N_terminal, C_terminal),
                    "hydro_moy": hydrophilie_moyenne(seq_peptide, proprietes_aa),
                    "coefficient_extinction": coefficient_extinction(seq_peptide),
                    "solubilite": solubilite(seq_peptide, N_terminal, C_terminal)

                }
                
                # aa, hydrophilicte, couleur = aa_vs_hydrophilicte(seq_peptide, proprietes_aa)
                # positions = list(range(1, len(seq_peptide) + 1))
                # fig_hydro = go.Figure()
                # fig_hydro.add_trace(go.Bar(
                #     x=positions,
                #     y=hydrophilicte,
                #     marker_color=couleur,
                #     text=aa,
                #     hovertemplate="AA: %{text}<br>Value: %{y}<extra></extra>"
                # ))
                # fig_hydro.update_layout(
                #     yaxis=dict(zeroline=True),
                #     xaxis=dict(
                #         dtick=1,
                #         tickmode="array",
                #         tickvals=positions,
                #         ticktext=aa,
                #     ),
                #     margin=dict(l=8, r=8, t=8, b=8),
                #     template="simple_white",
                # )

                graph_1_html = graphique_charge_ph(seq_peptide,N_terminal,C_terminal)
                graph_2_html = graphique_hydrophilicite(seq_peptide)
                #pio.to_html(fig_hydro, full_html=False)
    return render_template(
        "index.html",
        stat=stat,
        graph1=graph_1_html,
        graph2=graph_2_html,
        error_message=error_message,
        seq_peptide=seq_peptide,
        acid_amine=sorted(proprietes_aa.keys()),
        N_terminal=N_terminal,
        C_terminal=C_terminal,
    )

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=True)
