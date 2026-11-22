import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def cargar_y_analizar_3d(ruta_csv):
    print("Cargando datos para el análisis 3D...")
    df = pd.read_csv(ruta_csv)
    
    # Extraemos t, X (Re) e Y (Im)
    t = df.iloc[:, 0].values
    x = df.iloc[:, 1].values
    y = df.iloc[:, 2].values
    
    # Calculamos la velocidad de la parte real
    dx = np.gradient(x)
    
    # Creamos una figura con dos subplots 3D
    fig = make_subplots(rows=1, cols=2,
                        specs=[[{'type': 'scene'}, {'type': 'scene'}]],
                        subplot_titles=("Espacio de Fases 3D (Re, Im, dRe/dt)",
                                        "Sección de Poincaré 3D (t, Re, dRe/dt)"))

    # --- GRÁFICO 1: Trayectoria 3D Completa ---
    # Usamos el color para representar el tiempo (t)
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=dx,
        mode='lines',
        line=dict(width=1.5, color=t, colorscale='Plasma'),
        name='Trayectoria'
    ), row=1, col=1)

    # --- GRÁFICO 2: Sección de Poincaré 3D ---
    # Buscamos los cruces por cero de Im(Z)
    cruces = (y[:-1] < 0) & (y[1:] >= 0)
    t_c = t[:-1][cruces]
    x_c = x[:-1][cruces]
    dx_c = dx[:-1][cruces]

    fig.add_trace(go.Scatter3d(
        x=t_c, y=x_c, z=dx_c,
        mode='markers',
        marker=dict(size=2.5, color=t_c, colorscale='Plasma', opacity=0.8),
        name='Poincaré (Cruces)'
    ), row=1, col=2)

    # Configuración estética
    fig.update_layout(
        title="Buscando la Estructura Oculta de ζ(1/2 + it) en 3D",
        template="plotly_dark",
        height=700,
        showlegend=False,
        scene=dict(
            xaxis_title='Re(ζ)',
            yaxis_title='Im(ζ)',
            zaxis_title='d(Re(ζ))/dt',
            camera=dict(eye=dict(x=1.5, y=1.5, z=1.5)) # Vista isométrica
        ),
        scene2=dict(
            xaxis_title='t',
            yaxis_title='Re(ζ) en el corte',
            zaxis_title='d(Re(ζ))/dt',
            camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
        )
    )

    fig.write_html("atractor_3d.html")
    print("¡Listo! Abre atractor_3d.html")

if __name__ == "__main__":
    cargar_y_analizar_3d("salida_pfi01_N5000/trayectoria.csv")
