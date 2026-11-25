import sys
import numpy as np
from momentos_zeta import Z_rs      # mismo directorio: analisis_avanzado/

# Polinomio probado del 4o momento (Ivic-Conrey), segun Odlyzko, "zeta.moments.pdf", ec. (2.4):
#   int_0^T |zeta(1/2+it)|^4 dt = T * P2(ln(T/2pi)),   coeficientes de x^4 ... x^0:
C = [0.050660, 0.496227, 0.937279, 1.35334, -0.040924]
P = np.poly1d(C)
dP = P.deriv()


def media_local(L):
    # densidad local del momento acumulado: d/dT [T P2(ln(T/2pi))] = P2(L) + P2'(L)
    return P(L) + dP(L)


def main(Tmax=1e6):
    Ts = [x for x in (1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8) if x <= Tmax]
    print("M_2 = media de |zeta|^4 en [T, T+W], W=min(20000,T), paso 0.05; prediccion = polinomio de Ivic-Conrey\n")
    print(f"{'T':>10s} {'M_2 medido':>11s} {'+-':>8s} {'polinomio':>10s} {'medido/pol':>11s}")
    for T in Ts:
        W = min(20000.0, T)
        t = np.arange(T, T + W, 0.05)
        f = np.abs(Z_rs(t)) ** 4
        blk = np.array([b.mean() for b in np.array_split(f, 20)])
        se = blk.std(ddof=1) / np.sqrt(20)
        pred = np.mean(media_local(np.log(t / (2 * np.pi))))
        print(f"{T:10.0f} {f.mean():11.2f} {se:8.2f} {pred:10.2f} {f.mean() / pred:11.3f}")


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 1e6)
