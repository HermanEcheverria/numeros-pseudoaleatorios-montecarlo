"""
Parte 1: Mersenne Twister MT19937, mas el inciso (d).

El inciso (d) pide comparar un algoritmo propio contra `random` y `secrets`.
Se elige MT19937 porque es justamente el motor de `random`: si la
implementacion es correcta, la comparacion mide la calidad del codigo, no la
del algoritmo.

Produce:
    figuras/hist_mersenne.png
    resultados/mt_verificacion.tex
    resultados/mt_pruebas.tex
    resultados/comparacion_librerias.tex

    py -m scripts.p1_mersenne
"""

import random
import secrets
import time

from analisis import pruebas, salida
from prngs import mersenne
from semillas import N, SEMILLA_REFERENCIA_MT, SEMILLAS

SEMILLA = SEMILLAS["mersenne"]
N_VELOCIDAD = 200_000

REFERENCIA_5489 = [3499211612, 581869302, 3890346734, 3586334585, 545404204,
                   4161255391, 3922919429, 949333985, 2715962298, 1323567403]

ETIQUETAS = {"chi2": "$\\chi^2$", "ks": "K--S",
             "autocorrelacion": "Autocorr.", "rachas": "Rachas"}


def verificacion():
    """Contraste contra el vector publicado de la especificacion."""
    print("\n--- Verificacion contra el vector de referencia (OEIS A221557) ---")
    obtenido = mersenne.estados(10, SEMILLA_REFERENCIA_MT)
    coincide = obtenido == REFERENCIA_5489
    print(f"  semilla {SEMILLA_REFERENCIA_MT}, primeras 10 palabras de 32 bits")
    print(f"  coincide con la especificacion: {coincide}")

    filas = [[str(i + 1), str(obtenido[i]), str(REFERENCIA_5489[i]),
              "Sí" if obtenido[i] == REFERENCIA_5489[i] else "\\textbf{No}"]
             for i in range(5)]
    ruta = salida.guardar_tabla(
        ["$i$", "Implementación", "Especificación", "Coincide"],
        filas, "mt_verificacion.tex", alineacion="rrrc")
    print(f"  -> {ruta.name}")
    return coincide


def inciso_a_b():
    u = mersenne.uniformes(N, SEMILLA)
    print(f"\n  -> {salida.histogramas([('Mersenne Twister', u, SEMILLA)], 'hist_mersenne.png').name}")

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
        filas, "mt_pruebas.tex", alineacion="lrrl")
    print(f"  -> {ruta.name}")


def cronometrar(funcion, n):
    """Mejor de 3 corridas: el minimo filtra el ruido del sistema operativo."""
    return min(_una_corrida(funcion, n) for _ in range(3))


def _una_corrida(funcion, n):
    inicio = time.perf_counter()
    funcion(n)
    return time.perf_counter() - inicio


def inciso_d():
    """Comparacion contra random y secrets: velocidad y las cuatro pruebas."""
    print(f"\n--- (d) Comparacion con random y secrets (n = {N_VELOCIDAD:,}) ---")

    aleatorio = random.Random(SEMILLA)
    sistema = secrets.SystemRandom()

    fuentes = [
        ("MT19937 (propio)", lambda n: mersenne.uniformes(n, SEMILLA), "Sí"),
        ("random (MT19937 de C)", lambda n: [aleatorio.random() for _ in range(n)], "Sí"),
        ("secrets (SystemRandom)", lambda n: [sistema.random() for _ in range(n)], "\\textbf{No}"),
    ]

    filas = []
    for nombre, generar, reproducible in fuentes:
        segundos = cronometrar(generar, N_VELOCIDAD)
        por_millon = segundos * 1e6 / N_VELOCIDAD
        u = generar(N)
        resultados = pruebas.bateria(u)
        ps = " / ".join(salida.formato_p(p) for _, p in resultados.values())
        print(f"    {nombre:<24} {por_millon:>7.2f} s/10^6   reproducible: "
              f"{reproducible.replace(chr(92)+'textbf{','').replace('}','')}")
        print(f"      p-valores (chi2 / KS / autocorr. / rachas): {ps}")
        filas.append([nombre, f"{por_millon:.2f}", reproducible] +
                     [salida.formato_p(p) for _, p in resultados.values()])

    ruta = salida.guardar_tabla(
        ["Fuente", "s / $10^6$", "Reproducible", "$\\chi^2$", "K--S",
         "Autocorr.", "Rachas"],
        filas, "comparacion_librerias.tex", alineacion="lrcrrrr")
    print(f"\n  -> {ruta.name}")


def main():
    if not verificacion():
        raise SystemExit("La implementacion NO reproduce la especificacion.")
    inciso_a_b()
    inciso_d()
    print()


if __name__ == "__main__":
    main()
