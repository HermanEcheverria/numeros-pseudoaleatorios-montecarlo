"""
Parte 1, incisos (a) y (b) para LCG y RANDU, mas el estudio de robustez.

Produce:
    figuras/hist_lcg_randu.png     histogramas normalizados      (inciso a)
    figuras/pvalores.png           p-valores sobre K semillas    (segundo nivel)
    resultados/pruebas.tex         las 4 pruebas                 (inciso b)
    resultados/robustez.tex        tasas de rechazo y uniformidad

    py -m scripts.p1_pruebas
"""

from analisis import pruebas, robustez, salida
from prngs import lcg
from semillas import N, SEMILLAS

K = 100          # numero de semillas del estudio de robustez
BINS_CHI2 = 50   # clases para la prueba chi-cuadrado de uniformidad

GENERADORES = {
    "LCG": (lcg.PARAMS_NUMERICAL_RECIPES, SEMILLAS["lcg"]),
    "RANDU": (lcg.PARAMS_RANDU, SEMILLAS["randu"]),
}

ETIQUETAS = {
    "chi2": "$\\chi^2$",
    "ks": "K--S",
    "autocorrelacion": "Autocorr.",
    "rachas": "Rachas",
}


def secuencias():
    return {nombre: (lcg.uniformes(N, semilla, validar=False, **params), semilla)
            for nombre, (params, semilla) in GENERADORES.items()}


def inciso_a(seqs):
    print("\n--- (a) Histogramas ---")
    ruta = salida.histogramas(
        [(nombre, u, semilla) for nombre, (u, semilla) in seqs.items()],
        "hist_lcg_randu.png")
    print(f"  {ruta.name}")


def inciso_b(seqs):
    print("\n--- (b) Pruebas de hipotesis ---")
    filas = []
    for nombre, (u, semilla) in seqs.items():
        resultados = pruebas.bateria(u, k=BINS_CHI2)
        print(f"\n  {nombre} (semilla {semilla}, n = {N:,})")
        for prueba, (estadistico, p) in resultados.items():
            print(f"    {prueba:<16} estadistico = {estadistico:>10.4f}   "
                  f"p = {salida.formato_p(p):>9}   {pruebas.decision(p)}")
            filas.append([nombre, ETIQUETAS[prueba], f"{estadistico:.4f}",
                          salida.formato_p(p),
                          "No se rechaza" if p >= pruebas.ALPHA else "Se rechaza"])

    ruta = salida.guardar_tabla(
        ["Generador", "Prueba", "Estadístico", "$p$-valor", "Decisión ($\\alpha=0.05$)"],
        filas, "pruebas.tex", alineacion="llrrl")
    print(f"\n  -> {ruta.name}")


def estudio_robustez():
    print(f"\n--- Robustez: {K} semillas disjuntas por generador ---")
    estudios, filas = {}, []

    for nombre, (params, semilla) in GENERADORES.items():
        sems = robustez.semillas_disjuntas(semilla, K, N, **params)
        print(f"\n  {nombre}: semillas {sems[0]}, {sems[1]}, ... ({K} en total, "
              f"espaciadas {N:,} posiciones)")

        estudio = robustez.estudio(
            lambda n, s, p=params: lcg.uniformes(n, s, validar=False, **p),
            sems, N, k_bins=BINS_CHI2)
        estudios[nombre] = estudio

        for prueba, pvals in estudio.items():
            tasa, lo, hi, dentro = robustez.tasa_rechazo(pvals)
            chi2_u, p_u = robustez.uniformidad_pvalores(pvals)
            print(f"    {prueba:<16} rechazos = {tasa:5.1%} "
                  f"(banda NIST {lo:.1%}-{hi:.1%}) {'OK' if dentro else 'FUERA'}   "
                  f"uniformidad p-val: chi2 = {chi2_u:6.2f}, p = {salida.formato_p(p_u)}")
            filas.append([nombre, ETIQUETAS[prueba], f"{tasa:.1%}",
                          f"[{lo:.1%}, {hi:.1%}]", "Sí" if dentro else "\\textbf{No}",
                          f"{chi2_u:.2f}", salida.formato_p(p_u)])

    ruta = salida.guardar_tabla(
        ["Generador", "Prueba", "Rechazos", "Banda NIST", "Dentro",
         "$\\chi^2$ unif.", "$p$"],
        filas, "robustez.tex", alineacion="llrrcrr")
    print(f"\n  -> {ruta.name}")

    ruta = salida.histogramas_pvalores(estudios, "pvalores.png")
    print(f"  -> {ruta.name}")


def main():
    seqs = secuencias()
    inciso_a(seqs)
    inciso_b(seqs)
    estudio_robustez()
    print()


if __name__ == "__main__":
    main()
