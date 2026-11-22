import numpy as np
import plotly.graph_objects as go

def generar_agujero_negro(ruta_ceros):
    print("Cargando los 100,000 ceros de Odlyzko...")
    ceros = np.loadtxt(ruta_ceros)
    if ceros.ndim > 1:
        ceros = ceros[:, -1]
        
    print("Calculando espaciados...")
    # Calculamos las diferencias entre ceros consecutivos
    delta_n = np.diff(ceros)
    
    # Para ver el atractor, graficamos el espaciado actual vs el siguiente
    # x = Delta_n, y = Delta_{n+1}
    x = delta_n[:-1]
    y = delta_n[1:]
    z = np.arange(len(x)) # El eje Z será el índice (el paso del tiempo)
    
    print("Generando el mapa de retorno 3D...")
    
    fig = go.Figure(data=[go.Scatter3d(
        x=x, 
        y=y, 
        z=z,
        mode='markers',
        marker=dict(
            size=1.5, 
            color=z, # Color según el avance del tiempo
            colorscale='Twilight', # Una escala de colores preciosa
            opacity=0.5
        )
    )])

    # Ajustamos la cámara para mirar de frente al "agujero"
    fig.update_layout(
        title="Mapa de Retorno: El Atractor del Caos de Riemann (100k ceros)",
        scene=dict(
            xaxis_title='Espaciado actual (Δ_n)',
            yaxis_title='Espaciado siguiente (Δ_n+1)',
            zaxis_title='Tiempo (n)',
            camera=dict(
                eye=dict(x=0, y=0, z=2.5) # Vista cenital para ver la estructura 2D
            )
        ),
        template="plotly_dark",
        margin=dict(l=0, r=0, b=0, t=40)
    )
    
    fig.write_html("agujero_negro_real.html")
    print("¡Listo! Abre agujero_negro_real.html")

if __name__ == "__main__":
    generar_agujero_negro("ceros_100k.txt")
