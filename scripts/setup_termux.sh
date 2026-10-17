#!/data/data/com.termux/files/usr/bin/bash
set -e
echo "== Actualizando paquetes =="
pkg update -y && pkg upgrade -y
echo "== Instalando Python, git y numpy (paquete Termux) =="
pkg install -y python git python-numpy
echo "== Instalando dependencias de Python =="
pip install mpmath plotly
echo "== Listo. Proba con: python3 pfi01.py =="
