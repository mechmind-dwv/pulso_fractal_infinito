import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def cargar_y_analizar(ruta_csv):
    print("Cargando datos...")
    df = pd.read_csv(ruta_csv)
    
    # TRUCO DE INGENIERO: Extraer por posición, no por nombre.
    # Asumimos: Col 0 = t, Col 1 = Re(Z), Col 2 = Im(Z)
    # Esto evita el KeyError si las columnas tienen otros nombres.
    t = df.iloc[:, 0].values
    x = df.iloc[:, 1].values
    y = df.iloc[:, 2].values
    
    print(f"Datos cargados: {len(t)} muestras.")
    
    # 1. Derivadas numéricas (Velocidad)
    # np.gradient calcula la diferencia entre puntos consecutivos
    dx = np.gradient(x)
    
    # 2. Reconstrucción del espacio de fases (X vs dX)
    fig = make_subplots(rows=1, cols=2, 
                        subplot_titles=("Espacio de Fases (Re(Z) vs d(Re(Z)))", "Sección de Poincaré (Cruce de Im(Z) por 0)"))

    # Gráfico 1: Espacio de fases (Atractor)
    fig.add_trace(go.Scattergl(x=x, y=dx, mode='lines',
                               line=dict(width=0.5, color='cyan'),
                               name='Fase'), row=1, col=1)
    
    # Gráfico 2: Sección de Poincaré
    # Buscamos donde y cruza de negativo a positivo
    cruces = (y[:-1] < 0) & (y[1:] >= 0)
    t_cruces = t[:-1][cruces]
    x_cruces = x[:-1][cruces]
    
    fig.add_trace(go.Scattergl(x=t_cruces, y=x_cruces, mode='markers',
                               marker=dict(size=4, color='magenta', opacity=0.7),
                               name='Poincaré'), row=1, col=2)
    
    # Estética
    fig.update_layout(title="Buscando el Gran Atractor en la trayectoria de ζ(1/2 + it)",
                      template="plotly_dark",
                      showlegend=False,
                      height=600)
    fig.update_xaxes(title_text="Re(ζ)", row=1, col=1)
    fig.update_yaxes(title_text="d(Re(ζ))/dt", row=1, col=1)
    fig.update_xaxes(title_text="t", row=1, col=2)
    fig.update_yaxes(title_text="Re(ζ) en el corte", row=1, col=2)
    
    fig.write_html("atractor_zeta.html")
    print("¡Listo! Abre atractor_zeta.html para ver el atractor.")

if __name__ == "__main__":
    cargar_y_analizar("salida_pfi01_N5000/trayectoria.csv")
