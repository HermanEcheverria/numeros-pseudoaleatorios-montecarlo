"""
Parte 2, seccion 3.4: reduccion de varianza.

Se aplican variables antiteticas y variables de control a las integrales (a) y
(b) de la seccion 3.2. Las dos integrales estan elegidas a proposito: sin(pi x)
es SIMETRICA respecto a x = 1/2 y phi(x) es MONOTONA en [0,2], asi que cada
tecnica funciona en una y falla en la otra.

Produce:
    resultados/mc_varianza.tex

    py -m scripts.p2_varianza
"""

import math

import numpy as np

from analisis import salida
from montecarlo import varianza as vr
from prngs import mersenne
from semillas import SEMILLAS

N = 200_000


def f_a(x):
    return np.sin(np.pi * x)


def f_b(x):
    return np.exp(-x**2 / 2) / np.sqrt(2 * np.pi)


CASOS = {
    "(a) $\\sin(\\pi x)$ en $[0,1]$": (f_a, 0.0, 1.0, 2 / math.pi),
    "(b) $\\phi(x)$ en $[0,2]$": (f_b, 0.0, 2.0, 0.5 * math.erf(2 / math.sqrt(2))),
}

# Controles: g(x) = x es lineal;  g(x) = x(1-x) es simetrica como sin(pi x).
CONTROLES = {
    "$g(x)=x$": (lambda x: x, lambda a, b: (a + b) / 2),
    "$g(x)=x(1-x)$": (lambda x: x * (1 - x),
                      lambda a, b: (a + b) / 2 - (b**3 - a**3) / (3 * (b - a))),
}


def main():
    filas = []
    for etiqueta, (f, a, b, verdadero) in CASOS.items():
        limpio = etiqueta.split(" ")[0]
        print(f"\n{'=' * 72}\n{etiqueta}   valor real = {verdadero:.8f}\n{'=' * 72}")

        # --- Monte Carlo simple, la referencia ---------------------------
        u = mersenne.uniformes(N, SEMILLAS["mersenne"])
        muestra = vr.simple(f, a, b, u)
        I0, var0, se0, lo0, hi0 = vr.resumen(muestra)
        print(f"  {'simple':<22} I = {I0:.8f}  Var = {var0:.6e}  "
              f"EE = {se0:.2e}  factor = 1.00x")
        filas.append([limpio, "Simple", f"{I0:.6f}", f"{var0:.4e}",
                      "1.00", "1.00", f"[{lo0:.6f}, {hi0:.6f}]"])

        # --- Antiteticas: N/2 pares = N evaluaciones ---------------------
        u_mitad = mersenne.uniformes(N // 2, SEMILLAS["mersenne"])
        pares = vr.antiteticas(f, a, b, u_mitad)
        I1, var1, se1, lo1, hi1 = vr.resumen(pares)
        # Comparacion justa a igual numero de EVALUACIONES de f:
        # el estimador antitetico promedia N/2 pares, asi que su EE es
        # sqrt(var1/(N/2)); el factor es la razon de varianzas del promedio.
        factor1 = (var0 / N) / (var1 / (N // 2))
        print(f"  {'antitéticas':<22} I = {I1:.8f}  Var = {var1:.6e}  "
              f"EE = {se1:.2e}  factor = {factor1:.2f}x")
        filas.append([limpio, "Antitéticas", f"{I1:.6f}", f"{var1:.4e}",
                      f"{factor1:.2f}", f"{factor1:.2f}",
                      f"[{lo1:.6f}, {hi1:.6f}]"])

        # --- Variables de control ----------------------------------------
        for nombre, (g, media) in CONTROLES.items():
            media_g = media(a, b)
            corregida, c, rho = vr.control(f, a, b, u, g, media_g)
            I2, var2, se2, lo2, hi2 = vr.resumen(corregida)
            factor2 = var0 / var2
            print(f"  {'control ' + nombre:<22} I = {I2:.8f}  Var = {var2:.6e}  "
                  f"EE = {se2:.2e}  factor = {factor2:.2f}x   "
                  f"rho = {rho:+.4f}  c* = {c:+.4f}")
            filas.append([limpio, f"Control {nombre}", f"{I2:.6f}",
                          f"{var2:.4e}", f"{factor2:.2f}", f"{factor2:.2f}",
                          f"[{lo2:.6f}, {hi2:.6f}]"])

    ruta = salida.guardar_tabla(
        ["Integral", "Técnica", "$\\hat I$", "Varianza", "Factor",
         "\\emph{Speedup}", "IC 95\\%"],
        filas, "mc_varianza.tex", alineacion="llrrrrc")
    print(f"\n  -> {ruta.name}\n")


if __name__ == "__main__":
    main()
