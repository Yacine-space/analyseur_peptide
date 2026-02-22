import matplotlib.pyplot as plt
import numpy as np
from flask import Flask, request, render_template
import plotly.graph_objs as go
import plotly.io as pio


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

def ph(seq_peptide, proprietes_aa):
    ph = -1
    d={}
    
    while ph<=13:
        ph+=1
        for aa in seq_peptide:
            charge=infos_peptid(seq_peptide, proprietes_aa, ph)
            d[ph]=charge

    return d

def infos_peptid(seq_peptide, proprietes_aa,ph):

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
        
    return charge_nette

def masse_mol(seq_peptide, proprietes_aa):
    masse_mol_tot=0
    for aa in seq_peptide:
        if aa in proprietes_aa:
            masse_mol_tot+=proprietes_aa[aa]["masse_moleculaire"]

    masse_mol_pep=masse_mol_tot-(len(seq_peptide)-1)*18.015 #soustraire les liaisons hydro
    return masse_mol_pep

def calc_phi(seq_peptide, proprietes_aa):
    ph=0
    e={}
    while ph<14:
        ph+=0.01
        ph= round(ph, 2)
        phi=infos_peptid(seq_peptide, proprietes_aa, ph)
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

def charge_nette_ph_7(seq_peptide, proprietes_aa):
    charge_nette= ph(seq_peptide, proprietes_aa)
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

def charge_vs_ph(seq_peptide, proprietes_aa):
    valeurs_ph= np.arange(0, 14.1, 0.1)
    charges = [infos_peptid(seq_peptide, proprietes_aa, ph) for ph in valeurs_ph]
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
app = Flask(__name__)

@app.route('/', methods=["GET", "POST"])
def home():
    seq_peptide=""
    stat= None
    graph_1_html= ""
    graph_2_html = ""
    if request.method == "POST":
        seq_peptide= request.form.get("peptid", "").strip().upper()
        if seq_peptide:
            stat = {
                "letter_code_1": seq_peptide,
                'length': length_pep(seq_peptide),
                "letter_code_3": letter_code(seq_peptide, proprietes_aa),
                "poid_mol": masse_mol(seq_peptide, proprietes_aa),
                "phi": calc_phi(seq_peptide, proprietes_aa),
                "charge_nette_ph_7": charge_nette_ph_7(seq_peptide, proprietes_aa),
                "hydro_moy":hydrophilie_moyenne(seq_peptide, proprietes_aa)
            }
            valeurs_ph, charges = charge_vs_ph(seq_peptide, proprietes_aa)

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
                    dtick=1, #pas d e 1
                    tickmode="linear",
                    position= 0.5,
                    showline= True,
                    linecolor = "black",


                ),
                yaxis=dict(
                    title={
                        "text": "Charge nette",
                        "font": {"size": 14, "color": 'black'},
                    },
                    range=[-6, 6],
                    tickmode= "linear",
                    position= 0.5,
                    showline= True,
                    linecolor = "black",

                ),
                margin=dict(
                    l=8, 
                    r=8,    
                    t=8,   
                    b=8     
                ),
                width = 400,
                height=400    
       
            )


            aa, hydrophilicte, couleur = aa_vs_hydrophilicte(seq_peptide, proprietes_aa)
            positions = list(range(1, len(seq_peptide)))
            fig_hydro = go.Figure()
            fig_hydro.add_trace(go.Bar(
                x=positions,
                y=hydrophilicte,
                marker_color=couleur,
                text=aa,
                hovertemplate="AA: %{text}<br>Value: %{y}<extra></extra>"

                
                ))
            fig_hydro.update_layout(
                yaxis= dict(zeroline=True),
                xaxis=dict(
                    dtick=1,
                    tickmode="array",
                    tickvals = positions,
                    ticktext=aa,
                    ),
                margin=dict(
                    l=8, 
                    r=8,    
                    t=8,   
                    b=8     
                ),
                template="simple_white",

            )

            FloatingPointError
            graph_1_html = pio.to_html(fig_charge, full_html=False)
            graph_2_html = pio.to_html(fig_hydro, full_html=False)
    return(render_template("index.html", stat=stat, graph1=graph_1_html, graph2=graph_2_html))

if __name__=="__main__":
    app.run(debug=True)