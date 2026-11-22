import numpy as np
import plotly.graph_objects as go

def generar_red_cosmica_1M(ruta_ceros):
    print("Cargando el millón de ceros...")
    ceros = np.loadtxt(ruta_ceros)
    if ceros.ndim > 1:
        ceros = ceros[:, -1]
        
    print("Normalizando espaciados...")
    t = ceros[:-1]
    rho = (1 / (2 * np.pi)) * np.log(t / (2 * np.pi))
    delta_n = np.diff(ceros)
    s_n = delta_n * rho 
    
    # Muestreo de 100,000 puntos del millón para la gráfica
    indices = np.linspace(0, len(s_n) - 3, 100000, dtype=int)
    x = s_n[indices]
    y = s_n[indices + 1]
    z = s_n[indices + 2]
    
    densidad = np.abs(x - 1) + np.abs(y - 1) + np.abs(z - 1)
    
    print("Generando la Red Cósmica del Millón...")
    fig = go.Figure(data=[go.Scatter3d(
        x=x, y=y, z=z, mode='markers',
        marker=dict(size=1.2, color=densidad, colorscale='Electric', opacity=0.6)
    )])

    fig.update_layout(
        title="Red Cósmica de Riemann (Muestra de 100k del Millón)",
        scene=dict(
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            zaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            bgcolor='black', camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
        ),
        paper_bgcolor='black', plot_bgcolor='black', template="plotly_dark"
    )
    fig.write_html("red_cosmica_1M.html")
    print("¡Listo! Abre red_cosmica_1M.html")

if __name__ == "__main__":
    generar_red_cosmica_1M("ceros_1millon.txt")
