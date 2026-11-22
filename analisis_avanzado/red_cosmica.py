import numpy as np
import plotly.graph_objects as go

def generar_red_cosmica(ruta_ceros):
    print("Cargando los 100,000 ceros...")
    ceros = np.loadtxt(ruta_ceros)
    if ceros.ndim > 1:
        ceros = ceros[:, -1]
        
    print("Calculando la densidad y los espaciados normalizados...")
    # Normalización GUE
    t = ceros[:-1]
    rho = (1 / (2 * np.pi)) * np.log(t / (2 * np.pi))
    delta_n = np.diff(ceros)
    s_n = delta_n * rho 
    
    # Inmersión en 3D (Takens)
    # Para no colapsar tu móvil, tomamos una muestra de 40,000 puntos
    limite = 40000
    x = s_n[:limite]
    y = s_n[1:limite+1]
    z = s_n[2:limite+2]
    
    # Calculamos una "densidad local" para iluminar los filamentos
    # (Los sitios donde los espaciados son similares serán más brillantes)
    densidad = np.abs(x - 1) + np.abs(y - 1) + np.abs(z - 1)
    
    print("Generando la Red Cósmica (puede tardar unos segundos)...")
    
    fig = go.Figure(data=[go.Scatter3d(
        x=x, 
        y=y, 
        z=z,
        mode='markers',
        marker=dict(
            size=1.5, 
            color=densidad, # El color depende de la densidad local
            colorscale='Electric', # Escala de colores tipo nebulosa
            opacity=0.7,
            colorbar=dict(title="Densidad")
        )
    )])

    # Vista cenital para apreciar la red
    fig.update_layout(
        title="Red Cósmica de Riemann: La Estructura del Vacío",
        scene=dict(
            xaxis=dict(title='s_n', showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(title='s_{n+1}', showgrid=False, zeroline=False, showticklabels=False),
            zaxis=dict(title='s_{n+2}', showgrid=False, zeroline=False, showticklabels=False),
            bgcolor='black',
            camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
        ),
        paper_bgcolor='black',
        plot_bgcolor='black',
        template="plotly_dark",
        margin=dict(l=0, r=0, b=0, t=40)
    )
    
    fig.write_html("red_cosmica.html")
    print("¡Listo! Abre red_cosmica.html")

if __name__ == "__main__":
    generar_red_cosmica("ceros_100k.txt")
