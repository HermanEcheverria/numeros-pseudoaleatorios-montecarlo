"""
Parte 1: metodo de los cuadrados medios (von Neumann).

Produce:
    figuras/hist_cuadrados_medios.png
    resultados/cm_censo.tex      censo exhaustivo de las semillas de 4 digitos
    resultados/cm_pruebas.tex    las 4 pruebas de hipotesis

    py -m scripts.p1_cuadrados_medios
"""

from collections import Counter

from analisis import pruebas, salida
from prngs import middle_square as ms
from semillas import DIGITOS_CUADRADOS_MEDIOS as DIGITOS
from semillas import N, SEMILLAS

SEMILLA = SEMILLAS["cuadrados_medios"]


def verificacion():
    """Contraste contra dos secuencias publicadas, antes de usar el generador."""
    print("\n--- Verificacion contra secuencias publicadas ---")
    casos = [
        ("von Neumann, 2 dig., semilla 43", ms.estados(5, 43, 2), [84, 5, 2, 0, 0]),
        ("ciclo corto, 4 dig., semilla 6100", ms.estados(4, 6100, 4),
         [2100, 4100, 8100, 6100]),
    ]
    for etiqueta, obtenido, esperado in casos:
        print(f"  {etiqueta}: {obtenido} == {esperado} -> {obtenido == esperado}")


def censo_4_digitos():
    """Barrido exhaustivo de las 10,000 semillas de 4 digitos."""
    print("\n--- Censo exhaustivo, 4 digitos ---")
    orbitas = ms.censo(4)
    conteo = Counter(p for p, _ in orbitas.values())

    filas = [[str(p), f"{conteo[p]:,}", f"{100 * conteo[p] / 10_000:.1f}\\%"]
             for p in sorted(conteo)]
    for fila in filas:
        print(f"    periodo {fila[0]:>3}: {fila[1]:>6} semillas ({fila[2].replace(chr(92)+'%','%')})")
    print(f"    periodo maximo sobre las 10,000 semillas: {max(conteo)}")

    ruta = salida.guardar_tabla(
        ["Período", "Semillas", "Porcentaje"], filas, "cm_censo.tex",
        alineacion="rrr")
    print(f"  -> {ruta.name}")


def orbita_de_produccion():
    print(f"\n--- Semilla del informe: {SEMILLA} ({DIGITOS} digitos) ---")
    periodo, cola = ms.periodo_y_cola(SEMILLA, DIGITOS, limite=400_000)
    print(f"    cola    = {cola:,}   (estados que no se vuelven a visitar)")
    print(f"    periodo = {periodo:,}")
    print(f"    n = {N:,} < cola, asi que la muestra no contiene repeticiones")
    return periodo, cola


def inciso_a_b():
    u = ms.uniformes(N, SEMILLA, DIGITOS)

    ruta = salida.histogramas([("Cuadrados medios", u, SEMILLA)],
                             "hist_cuadrados_medios.png")
    print(f"\n  -> {ruta.name}")

    print("\n--- Pruebas de hipotesis ---")
    filas = []
    etiquetas = {"chi2": "$\\chi^2$", "ks": "K--S",
                 "autocorrelacion": "Autocorr.", "rachas": "Rachas"}
    for prueba, (estadistico, p) in pruebas.bateria(u).items():
        print(f"    {prueba:<16} estadistico = {estadistico:>10.4f}   "
              f"p = {salida.formato_p(p):>9}   {pruebas.decision(p)}")
        filas.append([etiquetas[prueba], f"{estadistico:.4f}",
                      salida.formato_p(p),
                      "No se rechaza" if p >= pruebas.ALPHA else "Se rechaza"])

    ruta = salida.guardar_tabla(
        ["Prueba", "Estadístico", "$p$-valor", "Decisión ($\\alpha=0.05$)"],
        filas, "cm_pruebas.tex", alineacion="lrrl")
    print(f"  -> {ruta.name}")


def main():
    verificacion()
    censo_4_digitos()
    orbita_de_produccion()
    inciso_a_b()
    print()


if __name__ == "__main__":
    main()
