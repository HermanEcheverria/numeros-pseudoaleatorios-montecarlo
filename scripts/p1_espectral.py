"""
Parte 1, inciso (c): prueba espectral visual para LCG y RANDU.

Produce:
    figuras/espectral_3d.png   ternas desde dos angulos, LCG vs RANDU
    figuras/espectral_k.png    histograma de k_i
    resultados/espectral.tex   planos contados

    py -m scripts.p1_espectral
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analisis import espectral, salida
from prngs import lcg
from semillas import SEMILLAS

N_PUNTOS = 3_000    # pocos puntos: con demasiados los planos se emborronan
N_CONTEO = 20_000   # muchos para contar planos y llenar el histograma

ANGULO_GENERICO = (30.0, -60.0)   # el que matplotlib usa por defecto


def secuencias(n):
    return {
        "LCG": lcg.uniformes(n, SEMILLAS["lcg"], **lcg.PARAMS_NUMERICAL_RECIPES),
        "RANDU": lcg.randu(n, SEMILLAS["randu"]),
    }


def figura_3d(seqs, angulo_canto):
    nombres = list(seqs)
    vistas = [("ángulo genérico", ANGULO_GENERICO),
              (f"de canto (elev {angulo_canto[0]:.0f}°, azim {angulo_canto[1]:.0f}°)",
               angulo_canto)]

    fig = plt.figure(figsize=(4.0 * len(nombres), 4.0 * len(vistas)))
    for i, (etiqueta_vista, (elev, azim)) in enumerate(vistas):
        for j, nombre in enumerate(nombres):
            eje = fig.add_subplot(len(vistas), len(nombres),
                                  i * len(nombres) + j + 1, projection="3d")
            x, y, z = espectral.ternas(seqs[nombre])
            eje.scatter(x, y, z, s=0.6, alpha=0.5, color="#4C72B0",
                        edgecolors="none")
            eje.view_init(elev=elev, azim=azim)
            eje.set_title(f"{nombre} — {etiqueta_vista}", fontsize=9)
            eje.set_xlabel("$u_i$", fontsize=7)
            eje.set_ylabel("$u_{i+1}$", fontsize=7)
            eje.set_zlabel("$u_{i+2}$", fontsize=7)
            eje.tick_params(labelsize=6)

    fig.suptitle(f"Ternas consecutivas, n = {N_PUNTOS:,} puntos", fontsize=10)
    fig.tight_layout()
    return salida._guardar(fig, "espectral_3d.png")


def figura_k(seqs):
    fig, ejes = plt.subplots(1, len(seqs), figsize=(4.4 * len(seqs), 3.0),
                             sharey=True)
    for eje, (nombre, u) in zip(ejes, seqs.items()):
        k = espectral.combinacion(u)
        eje.hist(k, bins=300, range=(-6, 10), color="#4C72B0")
        eje.set_title(f"{nombre}", fontsize=9)
        eje.set_xlabel("$k_i = 9u_i - 6u_{i+1} + u_{i+2}$")
        eje.set_xlim(-6, 10)
    ejes[0].set_ylabel("frecuencia")
    fig.suptitle("Distribución de $k_i$: RANDU solo alcanza 15 valores enteros",
                 fontsize=10)
    fig.tight_layout()
    return salida._guardar(fig, "espectral_k.png")


def main():
    (elev, azim), residuo = espectral.angulo_de_canto()
    print(f"\n--- (c) Prueba espectral ---")
    print(f"  Ángulo de canto derivado de d·h = 0: elev = {elev:.0f}°, "
          f"azim = {azim:.0f}°  (residuo {residuo:.1e})")

    seqs_conteo = secuencias(N_CONTEO)
    filas = []
    for nombre, u in seqs_conteo.items():
        k = espectral.combinacion(u)
        cuantos, valores, desviacion = espectral.contar_planos(k)
        print(f"\n  {nombre} (n = {N_CONTEO:,})")
        print(f"    planos con k entero  : {cuantos}")
        print(f"    valores de k         : {valores if cuantos else '—'}")
        print(f"    desviación máxima    : {desviacion:.2e}")
        filas.append([nombre, str(cuantos),
                      f"[{k.min():.3f}, {k.max():.3f}]", f"{desviacion:.1e}"])

    ruta = salida.guardar_tabla(
        ["Generador", "Planos ocupados", "Rango de $k_i$", "Desv. máx. al entero"],
        filas, "espectral.tex", alineacion="lrrr")
    print(f"\n  -> {ruta.name}")

    print(f"  -> {figura_3d(secuencias(N_PUNTOS), (elev, azim)).name}")
    print(f"  -> {figura_k(seqs_conteo).name}")
    print()


if __name__ == "__main__":
    main()
