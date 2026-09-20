"""
Parte 1: Blum Blum Shub.

Produce:
    figuras/hist_bbs.png
    resultados/bbs_condiciones.tex   condiciones sobre p y q
    resultados/bbs_pruebas.tex       las 4 pruebas de hipotesis

    py -m scripts.p1_bbs
"""

from analisis import pruebas, salida
from prngs import bbs
from semillas import N, SEMILLAS

SEMILLA = SEMILLAS["bbs"]

ETIQUETAS = {"chi2": "$\\chi^2$", "ks": "K--S",
             "autocorrelacion": "Autocorr.", "rachas": "Rachas"}


def condiciones():
    print("\n--- Condiciones sobre p y q ---")
    filas = []
    for descripcion, ok in bbs.verificar_primos():
        print(f"  [{'OK' if ok else 'NO'}] {descripcion}")
        filas.append([descripcion.replace("_", "\\_"),
                      "Sí" if ok else "\\textbf{No}"])

    ruta = salida.guardar_tabla(["Condición", "Se cumple"], filas,
                                "bbs_condiciones.tex", alineacion="lc")
    print(f"  -> {ruta.name}")


def parametros():
    periodo, lam = bbs.periodo_estados()
    print(f"\n--- Parametros ---")
    print(f"  M = p*q      = {bbs.M}  ({bbs.M.bit_length()} bits)")
    print(f"  lambda(M)    = {lam:,}")
    print(f"  periodo      = {periodo:,}  (~2^{periodo.bit_length() - 1})")
    print(f"  bits por paso= {bbs.bits_por_paso()} = floor(log2(log2(M)))")
    return periodo


def verificacion_periodo():
    """El periodo predicho debe coincidir con el observado en M pequeno."""
    print("\n--- Verificacion del periodo con primos seguros pequenos ---")
    for p, q, semilla in [(11, 23, 3), (23, 47, 5), (47, 59, 7)]:
        predicho, _ = bbs.periodo_estados(p, q)
        visto, observado = {}, None
        for i, x in enumerate(bbs.estados(400, semilla, p, q)):
            if x in visto:
                observado = i - visto[x]
                break
            visto[x] = i
        print(f"  p={p:>3} q={q:>3}: predicho {predicho:>4}, observado {observado:>4}"
              f"  -> {'coincide' if predicho == observado else 'NO COINCIDE'}")


def inciso_a_b():
    u = bbs.uniformes(N, SEMILLA)
    print(f"\n  -> {salida.histogramas([('Blum Blum Shub', u, SEMILLA)], 'hist_bbs.png').name}")

    print("\n--- Pruebas de hipotesis ---")
    filas = []
    for prueba, (estadistico, p) in pruebas.bateria(u).items():
        print(f"    {prueba:<16} estadistico = {estadistico:>10.4f}   "
              f"p = {salida.formato_p(p):>9}   {pruebas.decision(p)}")
        filas.append([ETIQUETAS[prueba], f"{estadistico:.4f}",
                      salida.formato_p(p),
                      "No se rechaza" if p >= pruebas.ALPHA else "Se rechaza"])

    ruta = salida.guardar_tabla(
        ["Prueba", "Estadístico", "$p$-valor", "Decisión ($\\alpha=0.05$)"],
        filas, "bbs_pruebas.tex", alineacion="lrrl")
    print(f"  -> {ruta.name}")


def main():
    condiciones()
    parametros()
    verificacion_periodo()
    inciso_a_b()
    print()


if __name__ == "__main__":
    main()
