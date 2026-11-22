import numpy as np

def calcular_li(ruta_ceros, max_n=100):
    print(f"Cargando ceros desde {ruta_ceros}...")
    ceros = np.loadtxt(ruta_ceros)
    if ceros.ndim > 1:
        ceros = ceros[:, -1]
        
    print(f"Calculando los primeros {max_n} coeficientes de Li...")
    lambdas = []
    
    # Fórmula de Keiper-Li usando los ceros complejos: rho = 1/2 + i*t
    for n in range(1, max_n + 1):
        # La suma es sobre todos los ceros y sus conjugados (por eso se multiplica por 2)
        # lambda_n = sum( 2 * (1 - Re( (1 - 1/rho)^n )) )
        # donde rho = 0.5 + i*t
        suma = 0.0
        for t in ceros:
            rho = 0.5 + 1j * t
            termino = 1 - (1 - 1/rho)**n
            suma += 2 * termino.real
            
        lambdas.append(suma)
        if n % 10 == 0:
            print(f"  Calculado n={n} -> lambda_{n} = {suma:.5f}")
            
    # Verificación
    print("\n" + "="*40)
    if all(l > 0 for l in lambdas):
        print("Todos los lambda_n > 0, pero es AUTOMATICO: con rho=1/2+it fijado cada termino 1-cos(n*phi) >= 0.")
        print("Esto NO prueba ni contrasta la Hipotesis (el resultado esta garantizado por construccion).")
    else:
        print("¡ALERTA! Se encontró un coeficiente negativo.")
        print("Esto sería un contraejemplo a la Hipótesis de Riemann.")
    print("="*40)
    
    # Guardar resultados
    np.savetxt("resultados_li.csv", np.column_stack((np.arange(1, max_n+1), lambdas)), 
               delimiter=",", header="n,lambda_n", comments="")
    print("Guardado en resultados_li.csv")

if __name__ == "__main__":
    # Usamos el archivo de 100k para que sea rápido
    calcular_li("ceros_100k.txt", max_n=100)
