import numpy as np
import plotly.graph_objects as go

def generar_atractor_masivo(ruta_ceros):
    print("Cargando datos masivos de ceros...")
    # Leemos el archivo de Odlyzko (asumiendo formato de una columna)
    try:
        ceros = np.loadtxt(ruta_ceros)
    except Exception as e:
        print(f"Error leyendo archivo: {e}")
        return

    # Si el archivo tiene varias columnas, tomamos la última (la altura t)
    if ceros.ndim > 1:
        ceros = ceros[:, -1]

    print(f"¡Cargados {len(ceros)} ceros!")
    
    # Nos quedamos con una muestra para que Plotly no colapse (ej. 50,000)
    # Plotly puede manejar 50k puntos en 3D si usamos Scatter3d y un buen móvil
    limite = min(50000, len(ceros) - 2)
    
    t_n = ceros[:limite]
    t_n1 = ceros[1:limite+1]
    delta_n = t_n1 - t_n
    
    print("Generando la bestia 3D...")
    
    fig = go.Figure(data=[go.Scatter3d(
        x=t_n, 
        y=t_n1, 
        z=delta_n,
        mode='markers',
        marker=dict(
            size=1, # Puntos muy pequeños
            color=delta_n, # El color depende del espaciado
            colorscale='Inferno', # Escala de colores tipo fuego/agujero negro
            opacity=0.6
        )
    )])

    fig.update_layout(
        title="Estructura del Agujero Negro de Riemann (100k ceros)",
        scene=dict(
            xaxis_title='t_n (Cero actual)',
            yaxis_title='t_{n+1} (Siguiente cero)',
            zaxis_title='Espaciado (Δ_n)'
        ),
        template="plotly_dark",
        margin=dict(l=0, r=0, b=0, t=40)
    )
    
    fig.write_html("atractor_agujero_negro.html")
    print("¡Listo! Abre atractor_agujero_negro.html")

if __name__ == "__main__":
    generar_atractor_masivo("ceros_100k.txt")
