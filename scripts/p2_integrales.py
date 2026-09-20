"""
Parte 2, seccion 3.2: integrales en una dimension.

    (a) int_0^1 sin(pi x) dx          = 2/pi
    (b) int_0^2 phi(x) dx             = Phi(2) - Phi(0)

Produce:
    figuras/convergencia.png       error absoluto vs N en log-log
    resultados/mc_integrales.tex   estimaciones con IC del 95%
    resultados/mc_pendientes.tex   pendientes ajustadas

    py -m scripts.p2_integrales
"""

import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analisis import robustez, salida
from montecarlo import estimador
from prngs import lcg, mersenne
from semillas import SEMILLAS

N_PRINCIPAL = 100_000
EXPONENTES = list(range(1, 7))     # N de 10^1 a 10^6
REPLICAS = 20                      # promediar el error para que la recta se vea


def f_a(x):
    return np.sin(np.pi * x)


def f_b(x):
    return np.exp(-x**2 / 2) / np.sqrt(2 * np.pi)


INTEGRALES = {
    "(a) $\\int_0^1 \\sin(\\pi x)\\,dx$": (f_a, 0.0, 1.0, 2 / math.pi),
    "(b) $\\int_0^2 \\phi(x)\\,dx$": (f_b, 0.0, 2.0,
                                     0.5 * math.erf(2 / math.sqrt(2))),
}

# Dos fuentes de uniformes: una libreria estandar y un PRNG propio (§3.2.3).
FUENTES = {
    "MT19937 (propio)": lambda n, s: mersenne.uniformes(n, s),
    "LCG (propio)": lambda n, s: lcg.uniformes(n, s, **lcg.PARAMS_NUMERICAL_RECIPES),
}


def estimaciones():
    print(f"\n--- Estimaciones con N = {N_PRINCIPAL:,} ---")
    filas = []
    for etiqueta, (f, a, b, verdadero) in INTEGRALES.items():
        for fuente, generar in FUENTES.items():
            u = generar(N_PRINCIPAL, SEMILLAS["mersenne"])
            I, se, lo, hi = estimador.integrar(f, a, b, u)
            atrapa = lo <= verdadero <= hi
            limpio = etiqueta.split(" ")[0]
            print(f"  {limpio} {fuente:<18} I = {I:.6f}  IC95 = [{lo:.6f}, {hi:.6f}]  "
                  f"real = {verdadero:.6f}  error = {I - verdadero:+.6f}  "
                  f"{'contiene' if atrapa else 'NO CONTIENE'}")
            filas.append([limpio, fuente, f"{I:.6f}", f"{se:.6f}",
                          f"[{lo:.6f}, {hi:.6f}]", f"{verdadero:.6f}",
                          f"{I - verdadero:+.2e}", "Sí" if atrapa else "\\textbf{No}"])

    ruta = salida.guardar_tabla(
        ["Integral", "Fuente", "$\\hat I_N$", "EE", "IC 95\\%", "Valor real",
         "Error", "Contiene"],
        filas, "mc_integrales.tex", alineacion="llrrcrrc")
    print(f"  -> {ruta.name}")


def _error_medio(f, a, b, verdadero, generar, n, replicas, semilla_base):
    """Error absoluto promedio sobre `replicas` bloques disjuntos de tamano n."""
    total = generar(n * replicas, semilla_base)
    errores = []
    for k in range(replicas):
        bloque = total[k * n:(k + 1) * n]
        I, *_ = estimador.integrar(f, a, b, bloque)
        errores.append(abs(I - verdadero))
    return float(np.mean(errores))


def convergencia():
    print(f"\n--- Estudio de convergencia ({REPLICAS} réplicas por N) ---")
    fig, ejes = plt.subplots(1, len(INTEGRALES), figsize=(5.2 * len(INTEGRALES), 3.8))
    filas = []

    for eje, (etiqueta, (f, a, b, verdadero)) in zip(ejes, INTEGRALES.items()):
        limpio = etiqueta.split(" ")[0]
        for fuente, generar in FUENTES.items():
            ns, errores = [], []
            for e in EXPONENTES:
                n = 10**e
                # con N grande se bajan las replicas para no reventar la memoria
                r = REPLICAS if e <= 4 else max(3, REPLICAS // (10 ** (e - 4)))
                errores.append(_error_medio(f, a, b, verdadero, generar, n, r,
                                            SEMILLAS["mersenne"]))
                ns.append(n)

            logn, loge = np.log10(ns), np.log10(errores)
            pendiente, intercepto = np.polyfit(logn, loge, 1)
            print(f"  {limpio} {fuente:<18} pendiente = {pendiente:+.4f}  "
                  f"(teórica -0.5, desvío {abs(pendiente + 0.5):.4f})")
            filas.append([limpio, fuente, f"{pendiente:+.4f}", "$-0.5$",
                          f"{abs(pendiente + 0.5):.4f}"])

            eje.loglog(ns, errores, "o-", markersize=4, label=f"{fuente}")

        # recta de referencia con pendiente exacta -1/2
        ref = np.array(ns, dtype=float)
        eje.loglog(ref, errores[0] * (ref / ref[0]) ** -0.5, "k--", linewidth=1,
                   label="$O(N^{-1/2})$")
        eje.set_title(etiqueta, fontsize=10)
        eje.set_xlabel("$N$")
        eje.set_ylabel("error absoluto medio")
        eje.legend(fontsize=7)
        eje.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    print(f"\n  -> {salida._guardar(fig, 'convergencia.png').name}")

    ruta = salida.guardar_tabla(
        ["Integral", "Fuente", "Pendiente ajustada", "Teórica", "Desvío"],
        filas, "mc_pendientes.tex", alineacion="llrrr")
    print(f"  -> {ruta.name}")


def verificacion_cobertura():
    """El IC del 95% debe atrapar el valor real ~95% de las veces."""
    print("\n--- Cobertura empírica del IC 95% ---")
    f, a, b, verdadero = INTEGRALES["(a) $\\int_0^1 \\sin(\\pi x)\\,dx$"]
    bloques, n = 1000, 10_000
    u = mersenne.uniformes(bloques * n, SEMILLAS["mersenne"])
    aciertos, _ = estimador.cobertura(f, a, b, u, bloques, verdadero)
    esperado = 0.95 * bloques
    ruido = math.sqrt(bloques * 0.95 * 0.05)
    print(f"  {aciertos}/{bloques} = {aciertos/bloques:.1%}  "
          f"(esperado {esperado:.0f} ± {ruido:.1f})")


def main():
    estimaciones()
    convergencia()
    verificacion_cobertura()
    print()


if __name__ == "__main__":
    main()
