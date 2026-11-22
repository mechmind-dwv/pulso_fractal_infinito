import numpy as np
import plotly.graph_objects as go

def analizar_100k(ruta_ceros):
    print("Cargando los 100,000 ceros de Odlyzko...")
    ceros = np.loadtxt(ruta_ceros)
    if ceros.ndim > 1:
        ceros = ceros[:, -1]
        
    print("Calculando espaciados normalizados...")
    # Espaciados
    t_n = ceros[:-1]
    delta_n = np.diff(ceros)
    
    # Densidad local (fórmula de Riemann)
    rho = (1 / (2 * np.pi)) * np.log(t_n / (2 * np.pi))
    
    # Espaciado normalizado (s_n)
    s_n = delta_n * rho
    
    # Filtramos valores extremos para que la gráfica no se vaya al infinito
    s_n = s_n[(s_n > 0) & (s_n < 4)]
    
    print(f"Analizando {len(s_n)} espaciados normalizados...")
    
    # 1. Calcular el histograma empírico (Densidad)
    hist_counts, bin_edges = np.histogram(s_n, bins=100, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # 2. Curvas teóricas
    s_teorico = np.linspace(0.01, 4, 200)
    # Wigner surmise para GUE: P(s) = (32/pi^2) * s^2 * exp(-4s^2/pi)
    p_gue = (32 / np.pi**2) * s_teorico**2 * np.exp(-4 * s_teorico**2 / np.pi)
    # Distribución de Poisson: P(s) = exp(-s)
    p_poisson = np.exp(-s_teorico)
    
    # 3. Estadístico <r> (Razón de espaciados consecutivos)
    # r_n = min(delta_n, delta_n+1) / max(delta_n, delta_n+1)
    d1 = delta_n[:-1]
    d2 = delta_n[1:]
    r_n = np.minimum(d1, d2) / np.maximum(d1, d2)
    r_promedio = np.mean(r_n)
    
    print("\n" + "="*50)
    print(f"ESTADÍSTICO <r> EMPÍRICO (100k ceros): {r_promedio:.4f}")
    print(f"Predicción GUE asintótica: 0.5996")
    print(f"Predicción Poisson asintótica: 0.386")
    print("="*50 + "\n")
    
    print("Generando la gráfica interactiva...")
    
    # 4. Construir la gráfica
    fig = go.Figure()
    
    # Histograma empírico (Barras o línea gruesa)
    fig.add_trace(go.Scatter(
        x=bin_centers, y=hist_counts,
        mode='lines',
        line=dict(color='#FFD700', width=3), # Dorado
        name=f'Empírico (100k ceros) | ⟨r⟩ = {r_promedio:.3f}'
    ))
    
    # Curva GUE
    fig.add_trace(go.Scatter(
        x=s_teorico, y=p_gue,
        mode='lines',
        line=dict(color='#00BFFF', width=3, dash='dash'), # Azul
        name='GUE (predicción cuántica)'
    ))
    
    # Curva Poisson
    fig.add_trace(go.Scatter(
        x=s_teorico, y=p_poisson,
        mode='lines',
        line=dict(color='#32CD32', width=3), # Verde
        name='Poisson (sin correlación)'
    ))
    
    fig.update_layout(
        title="N=100,000: Espaciados normalizados vs. GUE vs. Poisson",
        xaxis_title="s (espaciado / espaciado medio local)",
        yaxis_title="Densidad P(s)",
        template="plotly_dark",
        legend=dict(x=0.7, y=0.9, bgcolor="rgba(0,0,0,0.5)"),
        margin=dict(l=40, r=40, t=60, b=40)
    )
    
    fig.write_html("espectral_100k.html")
    print("¡Listo! Abre espectral_100k.html")

if __name__ == "__main__":
    analizar_100k("ceros_100k.txt")
