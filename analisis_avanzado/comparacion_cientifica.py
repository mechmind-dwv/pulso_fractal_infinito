import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def generar_comparacion(ruta_ceros):
    print("Cargando el millón de ceros (Datos empíricos)...")
    ceros = np.loadtxt(ruta_ceros)
    if ceros.ndim > 1:
        ceros = ceros[:, -1]
        
    # 1. Cálculo de espaciados normalizados empíricos
    t = ceros[:-1]
    rho = (1 / (2 * np.pi)) * np.log(t / (2 * np.pi))
    delta_n = np.diff(ceros)
    s_empirico = delta_n * rho 
    
    # Filtramos extremos para la gráfica 2D
    s_empirico_filtrado = s_empirico[(s_empirico > 0) & (s_empirico < 4)]
    
    # 2. Definiciones teóricas
    s_teorico = np.linspace(0.01, 4, 200)
    p_gue = (32 / np.pi**2) * s_teorico**2 * np.exp(-4 * s_teorico**2 / np.pi)
    p_poisson = np.exp(-s_teorico)
    
    # ==========================================
    # GRÁFICA 1: DISTRIBUCIÓN 2D (Cuantitativa)
    # ==========================================
    print("Generando gráfica 2D científica...")
    hist_counts, bin_edges = np.histogram(s_empirico_filtrado, bins=100, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    fig2d = go.Figure()
    fig2d.add_trace(go.Scatter(x=bin_centers, y=hist_counts, mode='lines',
                               line=dict(color='#FFD700', width=3), name='Empírico (Ceros de Riemann)'))
    fig2d.add_trace(go.Scatter(x=s_teorico, y=p_gue, mode='lines',
                               line=dict(color='#00BFFF', width=3, dash='dash'), name='Teórico GUE (Caos Cuántico)'))
    fig2d.add_trace(go.Scatter(x=s_teorico, y=p_poisson, mode='lines',
                               line=dict(color='#32CD32', width=3), name='Teórico Poisson (Azar)'))
    
    fig2d.update_layout(
        title="Distribución de Espaciados Normalizados P(s)",
        xaxis_title="s (espaciado / espaciado medio local)",
        yaxis_title="Densidad P(s)",
        template="plotly_dark", paper_bgcolor='black', plot_bgcolor='black',
        legend=dict(x=0.65, y=0.95, bgcolor="rgba(0,0,0,0.5)")
    )
    fig2d.write_html("cientifica_2d.html")

    # ==========================================
    # GRÁFICA 2: NUBES 3D (Geométrica)
    # ==========================================
    print("Generando gráfica 3D científica...")
    limite = 30000 # Muestra para la gráfica
    x_r = s_empirico[:limite]
    y_r = s_empirico[1:limite+1]
    z_r = s_empirico[2:limite+2]
    
    # Simulación Poisson
    s_poisson = np.random.exponential(scale=1.0, size=limite + 2)
    x_p = s_poisson[:limite]
    y_p = s_poisson[1:limite+1]
    z_p = s_poisson[2:limite+2]
    
    fig3d = make_subplots(rows=1, cols=2, 
                        specs=[[{'type': 'scene'}, {'type': 'scene'}]],
                        subplot_titles=("Empírico: Ceros de Riemann (GUE)", "Simulación: Azar Puro (Poisson)"))
    
    fig3d.add_trace(go.Scatter3d(x=x_r, y=y_r, z=z_r, mode='markers',
        marker=dict(size=1.5, color='#00BFFF', opacity=0.6), name='Riemann'), row=1, col=1)
        
    fig3d.add_trace(go.Scatter3d(x=x_p, y=y_p, z=z_p, mode='markers',
        marker=dict(size=1.5, color='#FF4500', opacity=0.6), name='Simulación'), row=1, col=2)

    # Ejes fijos en [0, 3] para una comparación justa
    axis_style = dict(range=[0, 3], showgrid=True, gridcolor='#333333', zeroline=False, showticklabels=False)
    
    fig3d.update_layout(
        title="El Duelo 3D: Estructura Cuántica vs. Ruido Aleatorio",
        scene=dict(xaxis=axis_style, yaxis=axis_style, zaxis=axis_style, bgcolor='black'),
        scene2=dict(xaxis=axis_style, yaxis=axis_style, zaxis=axis_style, bgcolor='black'),
        paper_bgcolor='black', plot_bgcolor='black', template="plotly_dark",
        showlegend=False
    )
    fig3d.write_html("cientifica_3d.html")
    print("¡Listo! Archivos generados: cientifica_2d.html y cientifica_3d.html")

if __name__ == "__main__":
    generar_comparacion("ceros_1millon.txt")
