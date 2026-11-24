#!/bin/bash
# Pipeline del ensayo 3. Ejecutar desde esta carpeta:  bash run_all.sh
# Datos: ../analisis_avanzado/ceros_100k.txt (primeros 100 000 ceros, t hasta ~74 920)
D=../analisis_avanzado/ceros_100k.txt
[ -f "$D" ] || { echo "Falta $D"; exit 1; }
{
  echo "=== Parte A: ceros de Z(t) con solo M enteros (altura 1e4) ==="
  python3 ceros_desde_enteros.py "$D" 10000 1500
  echo
  echo "=== Parte A: lo mismo a altura 6e4 ==="
  python3 ceros_desde_enteros.py "$D" 60000 1500
  echo
  echo "=== Parte B: picos de primos, ceros reales frente a GUE de igual densidad ==="
  python3 contraste_trazas.py "$D" 40000 400
} | tee resultados.txt
echo
echo "Salida guardada en resultados.txt"
