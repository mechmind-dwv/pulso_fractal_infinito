import numpy as np
import plotly.graph_objects as go

def generar_contraste(ruta_ceros):
    print("Cargando los 2 millones de ceros...")
    ceros = np.loadtxt(ruta_ceros)
    if ceros.ndim > 1:
        ceros = ceros[:, -1]
        
    # 1. Preparar datos de Riemann (GUE)
    t = ceros[:-1]
    rho = (1 / (2 * np.pi)) * np.log(t / (2 * np.pi))
    delta_n = np.diff(ceros)
    s_n = delta_n * rho 
    
    # Muestreo de 50k para la gráfica (Riemann)
    limite = 50000
    x_r = s_n[:limite]
    y_r = s_n[1:limite+1]
    z_r = s_n[2:limite+2]
    
    # 2. Generar datos Poisson (Azar puro)
    print("Generando el caos de Poisson (azar)...")
    # Los espaciados de Poisson siguen una distribución exponencial
    s_poisson = np.random.exponential(scale=1.0, size=limite + 2)
    x_p = s_poisson[:limite]
    y_p = s_poisson[1:limite+1]
    z_p = s_poisson[2:limite+2]
    
    print("Construyendo la gráfica de contraste...")
    
    fig = go.Figure()
    
    # Traza de Riemann (GUE - Estructura)
    fig.add_trace(go.Scatter3d(
        x=x_r, y=y_r, z=z_r, mode='markers',
        marker=dict(size=1.2, color='#00BFFF', opacity=0.8), # Azul
        name='Riemann (GUE) - Estructura'
    ))
    
    # Traza de Poisson (Caos - Ruido)
    fig.add_trace(go.Scatter3d(
        x=x_p, y=y_p, z=z_p, mode='markers',
        marker=dict(size=1.2, color='#FF4500', opacity=0.4), # Rojo/Naranja
        name='Poisson - Azar puro'
    ))

    fig.update_layout(
        title="GUE (Riemann) vs Poisson (Azar puro) en 3D",
        scene=dict(
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            zaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            bgcolor='black', camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
        ),
        paper_bgcolor='black', plot_bgcolor='black', template="plotly_dark",
        legend=dict(x=0.02, y=0.98)
    )
    
    fig.write_html("contraste_poisson.html")
    print("¡Listo! Abre contraste_poisson.html")

if __name__ == "__main__":
    generar_contraste("ceros_1millon.txt")
