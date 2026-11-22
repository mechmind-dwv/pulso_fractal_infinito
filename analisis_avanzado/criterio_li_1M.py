import numpy as np
import time

def calcular_li_1M(ruta_ceros, max_n=100):
    print(f"Cargando ceros desde {ruta_ceros}...")
    ceros = np.loadtxt(ruta_ceros)
    if ceros.ndim > 1:
        ceros = ceros[:, -1]
        
    print(f"Ceros cargados: {len(ceros)}. Calculando {max_n} coeficientes (Vectorizado)...")
    start = time.time()
    
    # Creamos los números complejos rho = 0.5 + i*t
    rho = 0.5 + 1j * ceros
    
    lambdas = []
    for n in range(1, max_n + 1):
        # np.power hace la potencia en paralelo para todos los ceros a la vez
        # (1 - 1/rho)^n
        terminos = 1 - np.power(1 - 1/rho, n)
        suma = np.sum(2 * terminos.real)
        lambdas.append(suma)
        
    end = time.time()
    print(f"¡Cálculo terminado en {end - start:.2f} segundos!")
    
    print("\n" + "="*40)
    print(f"Resultado para n={max_n}: lambda_{max_n} = {lambdas[-1]:.5f}")
    print(f"Resultado para n=1: lambda_1 = {lambdas[0]:.5f}")
    print("Todos positivos." if all(l > 0 for l in lambdas) else "¡ALERTA! Negativo encontrado.")
    print("="*40)

if __name__ == "__main__":
    calcular_li_1M("ceros_1millon.txt", max_n=100)
