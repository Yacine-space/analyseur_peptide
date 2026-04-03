from flask import Flask, request, render_template
from models.peptide import Peptide
from visualisation import plot
from config import AA_DATABASE
import os

app = Flask(__name__)

@app.route('/', methods=["GET", "POST"])
def home():
    stat= None
    error_message = None
    graph1=""
    graph2=""
    sequence = ""
    n_terminal=""
    c_terminal=""
    if request.method == "POST":
        sequence= request.form.get("peptid", "").strip().upper()
        n_terminal= request.form.get("N_term", "")
        c_terminal=request.form.get("C_term", "")

        try:
            p = Peptide(sequence, n_terminal, c_terminal)
            stat = {
                "letter_code_1": sequence,
                'length': p.longueur,
                "letter_code_3": p.letter_code,
                "poid_mol": p.masse_molaire,
                "phi": p.phi,
                "charge_nette_ph_7": p.charge_nette_ph_7,
                "hydro_moy": p.hydrophilie_moyenne,
                "coefficient_extinction": p.coefficient_extinction,
                "solubilite": p.solubilite

            }
                
            
            graph1=plot.graphique_charge_ph(p)
            graph2=plot.graphique_hydrophilicite(p)
        except ValueError as e:
            error_message = str(e)
    return render_template(
        "index.html",
        stat=stat,
        error_message=error_message,
        sequence=sequence,
        acid_amine=sorted(AA_DATABASE.keys()),
        n_terminal=n_terminal,
        c_terminal=c_terminal,
        graph1=graph1,
        graph2=graph2,
    )

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=True)
