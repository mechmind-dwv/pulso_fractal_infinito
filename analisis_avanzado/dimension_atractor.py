import numpy as np
import pandas as pd

def calcular_dimension_correlacion(ruta_csv):
    print("Calculando Dimensión de Correlación (D2)...")
    df = pd.read_csv(ruta_csv)
    
    # Extraemos los datos y calculamos la derivada
    t = df.iloc[:, 0].values
    x = df.iloc[:, 1].values
    y = df.iloc[:, 2].values
    dx = np.gradient(x)
    
    # 1. Extraemos la Sección de Poincaré
    cruces = (y[:-1] < 0) & (y[1:] >= 0)
    x_c = x[:-1][cruces]
    dx_c = dx[:-1][cruces]
    
    # Para el cálculo, usamos solo un subconjunto para que Termux no explote
    # Tomamos 2000 puntos distribuidos uniformemente
    indices = np.linspace(0, len(x_c)-1, min(2000, len(x_c)), dtype=int)
    puntos = np.column_stack((x_c[indices], dx_c[indices]))
    N = len(puntos)
    print(f"Analizando {N} puntos de la Sección de Poincaré...")
    
    # 2. Calculamos las distancias entre todos los pares de puntos
    distancias = []
    for i in range(N):
        for j in range(i+1, N):
            d = np.linalg.norm(puntos[i] - puntos[j])
            distancias.append(d)
            
    distancias = np.array(distancias)
    distancias = distancias[distancias > 0] # Quitamos distancias cero
    
    # 3. Función de Correlación C(r)
    # Contamos cuántos pares están a una distancia menor que 'r'
    radios = np.logspace(-2, 1, 20) # De 0.01 a 10
    C_r = []
    
    for r in radios:
        cuenta = np.sum(distancias < r) / len(distancias)
        C_r.append(cuenta)
        
    # 4. La pendiente de log(C(r)) vs log(r) es la dimensión fractal D2
    # Filtramos los valores válidos (entre 0 y 1)
    log_r = np.log(radios)
    log_C = np.log(C_r)
    
    # Regresión lineal en la zona de escalado (evitando los extremos)
    mascara = (log_C > -5) & (log_C < -0.1)
    pendiente, intercepto = np.polyfit(log_r[mascara], log_C[mascara], 1)
    
    print("\n" + "="*40)
    print(f"¡RESULTADO! Dimensión de Correlación D2 ≈ {pendiente:.4f}")
    print("="*40)
    if pendiente % 1 != 0:
        print("¡Es un número fractal! Firma de un Atractor Extraño.")
    else:
        print("Es un número entero. Podría ser ruido o una órbita regular.")
    print("(Para comparar: un plano 2D tiene D2=2. El atractor de Lorenz D2=2.06)")

if __name__ == "__main__":
    calcular_dimension_correlacion("salida_pfi01_N5000/trayectoria.csv")
