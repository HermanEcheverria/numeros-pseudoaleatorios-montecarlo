"""
Parte 2, seccion 3.3: alta dimension y la maldicion de la dimensionalidad.

Produce:
    figuras/fraccion_volumen.png     V_d / 2^d contra d
    resultados/mc_bola.tex           estimaciones con IC para d = 2,5,10,20
    resultados/mc_rejilla.tex        puntos que necesitaria la rejilla

    py -m scripts.p2_alta_dimension
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analisis import salida
from montecarlo import alta_dimension as ad
from prngs import mersenne
from semillas import SEMILLAS

DIMENSIONES = [2, 5, 10, 20]
N_PUNTOS = 200_000
M_REJILLA = 10


def estimaciones():
    print(f"\n--- Volumen de la bola unitaria, N = {N_PUNTOS:,} puntos ---")
    filas = []
    for d in DIMENSIONES:
        u = mersenne.uniformes(N_PUNTOS * d, SEMILLAS["mersenne"])
        V, se, lo, hi, N, aciertos = ad.monte_carlo(d, u)
        esperados = N * ad.volumen_exacto(d) / 2.0**d
        real = ad.volumen_exacto(d)
        atrapa = lo <= real <= hi
        error_rel = abs(V - real) / real if real else float("nan")
        print(f"  d = {d:>2}:  aciertos = {aciertos:>7,} (esperados "
              f"{esperados:>9.3f})  V = {V:.6f}  IC95 = [{lo:.6f}, {hi:.6f}]  "
              f"real = {real:.6f}  err.rel = {error_rel:.1%}")
        filas.append([str(d), f"{aciertos:,}".replace(",", "\\,"),
                      f"{esperados:.3f}", f"{V:.6f}",
                      f"[{lo:.6f}, {hi:.6f}]", f"{real:.6f}",
                      f"{error_rel:.2%}", "Sí" if atrapa else "\\textbf{No}"])

    ruta = salida.guardar_tabla(
        ["$d$", "Aciertos", "Esperados", "$\\hat V_d$", "IC 95\\%",
         "$V_d$ exacto", "Error rel.", "Contiene"],
        filas, "mc_bola.tex", alineacion="rrrrcrrc")
    print(f"  -> {ruta.name}")


def rejilla():
    print(f"\n--- Cuadratura en rejilla con m = {M_REJILLA} por eje ---")
    filas = []
    for d in DIMENSIONES:
        puntos = ad.puntos_para_rejilla(d, M_REJILLA)
        V_rejilla, _ = ad.cuadratura_rejilla(d, M_REJILLA)
        real = ad.volumen_exacto(d)

        if V_rejilla is None:
            veredicto = "inviable"
            err = "---"
            print(f"  d = {d:>2}: {puntos:>25,} evaluaciones  -> INVIABLE")
        else:
            err = f"{abs(V_rejilla - real) / real:.2%}"
            veredicto = f"{V_rejilla:.6f}"
            print(f"  d = {d:>2}: {puntos:>25,} evaluaciones  "
                  f"V = {V_rejilla:.6f}  err.rel = {err}")

        # comparacion a igual presupuesto: MC con el mismo numero de puntos
        filas.append([str(d), f"${M_REJILLA}^{{{d}}}$", f"{puntos:,}".replace(",", "\\,"),
                      veredicto, err])

    ruta = salida.guardar_tabla(
        ["$d$", "Puntos", "Evaluaciones", "$\\hat V_d$ rejilla", "Error rel."],
        filas, "mc_rejilla.tex", alineacion="rlrrr")
    print(f"  -> {ruta.name}")

    print(f"\n  Con el mismo presupuesto de {N_PUNTOS:,} puntos, Monte Carlo")
    print(f"  entrega una estimacion en TODA dimension; la rejilla ya no llega")
    print(f"  en d = 10, donde necesitaria 10^10 = 10,000 millones.")


def comparacion_a_igual_presupuesto():
    """Error de MC vs rejilla con el mismo numero de evaluaciones."""
    print("\n--- Mismo presupuesto: MC contra rejilla ---")
    filas = []
    for d in DIMENSIONES:
        # m tal que m^d sea lo mas cercano posible a N_PUNTOS sin pasarse
        m = max(2, int(N_PUNTOS ** (1.0 / d)))
        presupuesto = m**d
        V_rej, _ = ad.cuadratura_rejilla(d, m)
        u = mersenne.uniformes(presupuesto * d, SEMILLAS["mersenne"])
        V_mc, *_ = ad.monte_carlo(d, u)
        real = ad.volumen_exacto(d)

        err_rej = abs(V_rej - real) / real if V_rej is not None else None
        err_mc = abs(V_mc - real) / real
        texto_rej = f"{err_rej:.2%}" if err_rej is not None else "inviable"
        print(f"  d = {d:>2}: m = {m:>3}, presupuesto = {presupuesto:>9,}   "
              f"rejilla {texto_rej:>9}   MC {err_mc:>8.2%}")
        filas.append([str(d), str(m), f"{presupuesto:,}".replace(",", "\\,"),
                      texto_rej, f"{err_mc:.2%}"])

    ruta = salida.guardar_tabla(
        ["$d$", "$m$", "Presupuesto", "Error rejilla", "Error MC"],
        filas, "mc_presupuesto.tex", alineacion="rrrrr")
    print(f"  -> {ruta.name}")


def figura_fraccion():
    ds = np.arange(1, 31)
    fraccion = np.array([ad.volumen_exacto(d) / 2.0**d for d in ds])

    fig, eje = plt.subplots(figsize=(6.2, 3.8))
    eje.semilogy(ds, fraccion, "o-", markersize=4, color="#4C72B0")
    for d in DIMENSIONES:
        eje.annotate(f"$d={d}$\n{ad.volumen_exacto(d)/2.0**d:.2e}",
                     (d, ad.volumen_exacto(d) / 2.0**d),
                     textcoords="offset points", xytext=(6, 6), fontsize=7)
    eje.set_xlabel("dimensión $d$")
    eje.set_ylabel("$V_d / 2^d$  (escala log)")
    eje.set_title("Fracción del cubo que ocupa la bola unitaria")
    eje.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    print(f"\n  -> {salida._guardar(fig, 'fraccion_volumen.png').name}")

    print("\n  Fracción del cubo ocupada por la bola:")
    for d in DIMENSIONES + [30]:
        print(f"    d = {d:>2}: {ad.volumen_exacto(d) / 2.0**d:.3e}")


def main():
    estimaciones()
    rejilla()
    comparacion_a_igual_presupuesto()
    figura_fraccion()
    print()


if __name__ == "__main__":
    main()
