from config import AA_DATABASE
import plotly.graph_objs as go
import plotly.io as pio
from services import analyses

def graphique_charge_ph(peptide):
    valeurs_ph, charges = analyses.charge_vs_ph(peptide)
    fig_charge = go.Figure()
    fig_charge.add_trace(go.Scatter(
        x=valeurs_ph,
        y=charges,
        mode="lines",
        line=dict(color="darkblue", width=3),
    ))

    fig_charge.update_layout(
        xaxis=dict(
            range=[0, 14],
            dtick=1,
            tickmode="linear",
            position=0.5,   
            showline=True,
            linecolor="black",
        ),
        yaxis=dict(
            range=[-6, 6],
            tickmode="linear",
            position=0.5,
            showline=True,
            linecolor="black",
        ),
        margin=dict(l=10, r=10, t=10, b=10),
        width=400,
        height=400,
    )
    fig_charge.add_annotation(
        x=1, y=0.5,
        xref="paper", yref="paper",
        text="<b>pH</b>",
        showarrow=False,
        xanchor="left",
        yanchor="middle"
    )
    fig_charge.add_annotation(
        x=0.5, y=1.06,
        xref="paper", yref="paper",
        text="<b>Charge nette</b>",
        showarrow=False,
        xanchor="center",
        yanchor="top",
    )
    graph_1_html = pio.to_html(fig_charge, full_html=False)
    return graph_1_html

def graphique_hydrophilicite(peptide):
    aa, hydrophilicte, couleur = analyses.aa_vs_hydrophilicte(peptide)
    positions = list(range(1, len(peptide.sequence) + 1))
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
        template="simple_white",
    )
    graph_2_html = pio.to_html(fig_hydro, full_html=False)
    return graph_2_html