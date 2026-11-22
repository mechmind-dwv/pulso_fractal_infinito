import numpy as np
import plotly.graph_objects as go

def generar_atractor_gue(ruta_ceros):
    print("Cargando los 100,000 ceros...")
    ceros = np.loadtxt(ruta_ceros)
    if ceros.ndim > 1:
        ceros = ceros[:, -1]
        
    print("Calculando la densidad local y normalizando...")
    # La densidad local de ceros es aproximadamente (1/2pi) * ln(t/2pi)
    t = ceros[:-1]
    rho = (1 / (2 * np.pi)) * np.log(t / (2 * np.pi))
    
    # Espaciados reales y normalizados
    delta_n = np.diff(ceros)
    s_n = delta_n * rho # ¡Aquí está la magia! Espaciado normalizado
    
    # Mapa de retorno normalizado
    x = s_n[:-1]
    y = s_n[1:]
    z = np.arange(len(x))
    
    print("Generando el Atractor Universal...")
    
    fig = go.Figure(data=[go.Scatter3d(
        x=x, 
        y=y, 
        z=z,
        mode='markers',
        marker=dict(
            size=1.2, 
            color=z, 
            colorscale='Jet', 
            opacity=0.6
        )
    )])

    fig.update_layout(
        title="Mapa de Retorno Normalizado (GUE): La Estructura Universal",
        scene=dict(
            xaxis_title='Espaciado normalizado actual (s_n)',
            yaxis_title='Espaciado normalizado siguiente (s_n+1)',
            zaxis_title='Tiempo (n)',
            camera=dict(eye=dict(x=0, y=0, z=2.5))
        ),
        template="plotly_dark",
        margin=dict(l=0, r=0, b=0, t=40)
    )
    
    fig.write_html("atractor_gue.html")
    print("¡Listo! Abre atractor_gue.html")

if __name__ == "__main__":
    generar_atractor_gue("ceros_100k.txt")
