import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def generar_contraste_visual(ruta_ceros):
    print("Cargando los 2 millones de ceros...")
    ceros = np.loadtxt(ruta_ceros)
    if ceros.ndim > 1:
        ceros = ceros[:, -1]
        
    # 1. Preparar datos de Riemann (GUE)
    t = ceros[:-1]
    rho = (1 / (2 * np.pi)) * np.log(t / (2 * np.pi))
    delta_n = np.diff(ceros)
    s_n = delta_n * rho 
    
    limite = 40000 # 40k puntos para que la gráfica vuele
    x_r = s_n[:limite]
    y_r = s_n[1:limite+1]
    z_r = s_n[2:limite+2]
    
    # 2. Generar datos Poisson (Azar puro)
    print("Generando el caos de Poisson...")
    s_poisson = np.random.exponential(scale=1.0, size=limite + 2)
    x_p = s_poisson[:limite]
    y_p = s_poisson[1:limite+1]
    z_p = s_poisson[2:limite+2]
    
    print("Construyendo la gráfica comparativa...")
    
    # Creamos dos escenas 3D lado a lado
    fig = make_subplots(rows=1, cols=2, 
                        specs=[[{'type': 'scene'}, {'type': 'scene'}]],
                        subplot_titles=("Riemann (GUE) - Estructura", "Poisson - Azar Puro"))
    
    # Traza de Riemann (Izquierda)
    fig.add_trace(go.Scatter3d(
        x=x_r, y=y_r, z=z_r, mode='markers',
        marker=dict(size=1.5, color='#00BFFF', opacity=0.5),
        name='Riemann (GUE)'
    ), row=1, col=1)
    
    # Traza de Poisson (Derecha)
    fig.add_trace(go.Scatter3d(
        x=x_p, y=y_p, z=z_p, mode='markers',
        marker=dict(size=1.5, color='#FF4500', opacity=0.5),
        name='Poisson (Azar)'
    ), row=1, col=2)

    # FORZAMOS LOS LÍMITES DE LOS EJES A [0, 3] PARA AMBOS
    # Esto corta las colas largas de Poisson y nos deja ver el núcleo
    axis_style = dict(
        range=[0, 3], 
        showgrid=False, zeroline=False, showticklabels=False,
        backgroundcolor="black", showbackground=True
    )

    fig.update_layout(
        title="El Duelo: Orden Cuántico vs. Ruido Aleatorio",
        scene=dict(xaxis=axis_style, yaxis=axis_style, zaxis=axis_style, bgcolor='black'),
        scene2=dict(xaxis=axis_style, yaxis=axis_style, zaxis=axis_style, bgcolor='black'),
        paper_bgcolor='black', plot_bgcolor='black', template="plotly_dark",
        margin=dict(l=0, r=0, b=0, t=50)
    )
    
    fig.write_html("contraste_visual.html")
    print("¡Listo! Abre contraste_visual.html")

if __name__ == "__main__":
    generar_contraste_visual("ceros_1millon.txt")
