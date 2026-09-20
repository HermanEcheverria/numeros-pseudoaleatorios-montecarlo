"""
Parte 2, seccion 3.5: cuando el generador arruina el resultado.

Se estima por Monte Carlo la integral triple

    I = int_[0,1]^3 1[x^2+y^2+z^2 <= 1] dx dy dz = pi/6

usando ternas consecutivas, que es donde vive el defecto de RANDU, y se repite
con Mersenne Twister. Despues se muestra a que escala aparece realmente el
fallo.

Produce:
    figuras/randu_convergencia.png
    resultados/mc_randu.tex
    resultados/mc_randu_escala.tex

    py -m scripts.p2_randu
"""

import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analisis import salida
from prngs import lcg, mersenne
from semillas import SEMILLAS

VERDADERO = math.pi / 6
N_TERNAS = 2_000_000
N_ESCALA = 1_000_000
Z_95 = 1.959963984540054

# Espaciado entre los 15 planos de RANDU: la normal es h = (9,-6,1) y dos
# planos consecutivos distan 1/|h|.
ESPACIADO = 1.0 / math.sqrt(9**2 + 6**2 + 1**2)

FUENTES = {
    "RANDU": lambda n: lcg.randu(n, SEMILLAS["randu"]),
    "Mersenne Twister": lambda n: mersenne.uniformes(n, SEMILLAS["mersenne"]),
    "LCG (Num. Recipes)": lambda n: lcg.uniformes(
        n, SEMILLAS["lcg"], **lcg.PARAMS_NUMERICAL_RECIPES),
}


def _puntos(u, n):
    return np.asarray(u[:n * 3], dtype=float).reshape(n, 3)


def estimar(u, n_ternas):
    """Fraccion de ternas disjuntas dentro del octante de la bola unitaria."""
    dentro = np.sum(_puntos(u, n_ternas) ** 2, axis=1) <= 1.0
    p = float(dentro.mean())
    se = math.sqrt(p * (1 - p) / n_ternas)
    return p, se, p - Z_95 * se, p + Z_95 * se


def tabla():
    print(f"\n--- Integral triple, {N_TERNAS:,} ternas ---")
    print(f"  valor real = pi/6 = {VERDADERO:.8f}")
    filas = []
    for nombre, generar in FUENTES.items():
        I, se, lo, hi = estimar(generar(N_TERNAS * 3), N_TERNAS)
        error = I - VERDADERO
        sigmas = abs(error) / se
        atrapa = lo <= VERDADERO <= hi
        print(f"  {nombre:<20} I = {I:.8f}  error = {error:+.2e} = "
              f"{sigmas:4.1f} EE   {'contiene' if atrapa else 'NO CONTIENE'}")
        filas.append([nombre, f"{I:.6f}", f"[{lo:.6f}, {hi:.6f}]",
                      f"{error:+.2e}", f"{sigmas:.1f}",
                      "Sí" if atrapa else "\\textbf{No}"])

    ruta = salida.guardar_tabla(
        ["Fuente", "$\\hat I$", "IC 95\\%", "Error", "$|$error$|/$EE",
         "Contiene $\\pi/6$"],
        filas, "mc_randu.tex", alineacion="lrcrrc")
    print(f"  -> {ruta.name}")
    print("\n  Con una region grande y suave los tres aciertan: los 15 planos")
    print("  atraviesan la bola muchas veces y actuan como estratificacion.")


def escala():
    """
    El fallo de RANDU aparece cuando la region mide menos que el espaciado
    entre planos. Se estima el volumen de bolas centradas y cada vez mas
    pequenas.
    """
    print("\n--- Bolas pequenas: donde el defecto si muerde ---")
    print(f"  espaciado entre planos de RANDU = 1/|h| = {ESPACIADO:.4f}")

    fuentes = {n: g for n, g in FUENTES.items() if n != "LCG (Num. Recipes)"}
    muestras = {nombre: _puntos(generar(N_ESCALA * 3), N_ESCALA)
                for nombre, generar in fuentes.items()}

    cabecera = f"  {'radio':>7} {'V real':>11}"
    for nombre in muestras:
        cabecera += f" {nombre[:8]:>11} {'err':>8}"
    print(cabecera)

    filas = []
    for r in (0.30, 0.15, 0.08, 0.05, 0.03):
        real = 4 / 3 * math.pi * r**3
        linea = f"  {r:>7.2f} {real:>11.3e}"
        fila = [f"{r:.2f}", f"{real:.3e}"]
        for puntos in muestras.values():
            est = float((np.sum((puntos - 0.5) ** 2, axis=1) <= r * r).mean())
            err = abs(est - real) / real
            linea += f" {est:>11.3e} {err:>7.1%}"
            fila += [f"{est:.3e}", f"{err:.1%}"]
        print(linea)
        filas.append(fila)

    encabezados = ["$r$", "$V$ real"]
    for nombre in muestras:
        encabezados += [nombre, "Error"]
    ruta = salida.guardar_tabla(encabezados, filas, "mc_randu_escala.tex",
                                alineacion="rr" + "rr" * len(muestras))
    print(f"  -> {ruta.name}")


def convergencia():
    """Error contra N: con sesgo, la curva se estanca en vez de bajar."""
    tamanos = [10**e for e in range(2, 7)]
    fig, ejes = plt.subplots(1, 2, figsize=(11, 4.0))

    for nombre, generar in FUENTES.items():
        u = generar(max(tamanos) * 3)
        ejes[0].loglog(tamanos, [abs(estimar(u, n)[0] - VERDADERO) for n in tamanos],
                       "o-", markersize=4, label=nombre)

    ref = np.array(tamanos, dtype=float)
    ejes[0].loglog(ref, 0.2 * ref**-0.5, "k--", linewidth=1, label="$O(N^{-1/2})$")
    ejes[0].set_xlabel("número de ternas $N$")
    ejes[0].set_ylabel(r"$|\hat I_N - \pi/6|$")
    ejes[0].set_title("Bola unitaria: los tres convergen")
    ejes[0].legend(fontsize=7)
    ejes[0].grid(True, which="both", alpha=0.3)

    radios = np.array([0.30, 0.20, 0.15, 0.10, 0.08, 0.05, 0.03])
    for nombre, generar in FUENTES.items():
        if nombre == "LCG (Num. Recipes)":
            continue
        puntos = _puntos(generar(N_ESCALA * 3), N_ESCALA)
        errores = []
        for r in radios:
            real = 4 / 3 * math.pi * r**3
            est = float((np.sum((puntos - 0.5) ** 2, axis=1) <= r * r).mean())
            errores.append(abs(est - real) / real)
        ejes[1].loglog(radios, errores, "o-", markersize=4, label=nombre)

    ejes[1].axvline(ESPACIADO, color="#C44E52", linestyle=":", linewidth=1.5,
                    label=f"espaciado = {ESPACIADO:.3f}")
    ejes[1].set_xlabel("radio de la bola")
    ejes[1].set_ylabel("error relativo")
    ejes[1].set_title("Bolas pequeñas: RANDU se despeña bajo el espaciado")
    ejes[1].legend(fontsize=7)
    ejes[1].grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    print(f"\n  -> {salida._guardar(fig, 'randu_convergencia.png').name}")


def main():
    tabla()
    escala()
    convergencia()
    print()


if __name__ == "__main__":
    main()
