"""
Parte 1, requisito 2: tabla comparativa de los cinco generadores.

Velocidad y periodo NO se copian de ninguna fuente: la velocidad se cronometra
en esta maquina y el periodo sale del modulo de cada generador (derivado donde
hay teoria, medido donde no la hay, que es el caso de cuadrados medios).

Produce:
    resultados/tabla_comparativa.tex

    py -m scripts.p1_tabla_comparativa
"""

import time

from analisis import salida
from prngs import bbs, lcg, mersenne
from prngs import middle_square as ms
from semillas import DIGITOS_CUADRADOS_MEDIOS as DIGITOS
from semillas import SEMILLAS

N_VELOCIDAD = 50_000


def cronometrar(generar, n=N_VELOCIDAD):
    """Mejor de 3: el minimo filtra el ruido del sistema operativo."""
    tiempos = []
    for _ in range(3):
        inicio = time.perf_counter()
        generar(n)
        tiempos.append(time.perf_counter() - inicio)
    return min(tiempos) * 1e6 / n   # segundos por millon


GENERADORES = [
    {
        "nombre": "LCG (Num. Recipes)",
        "generar": lambda n: lcg.uniformes(n, SEMILLAS["lcg"],
                                           **lcg.PARAMS_NUMERICAL_RECIPES),
        "periodo": "$2^{32}$",
        "origen": "Hull--Dobell",
        "seguridad": "Nula: la forma cerrada permite recuperar $(a,c,m)$ "
                     "con pocas salidas.",
        "usos": "Simulación sencilla, videojuegos, sistemas embebidos, "
                "\\texttt{rand()} de C.",
    },
    {
        "nombre": "RANDU",
        "generar": lambda n: lcg.randu(n, SEMILLAS["randu"]),
        "periodo": "$2^{29}$",
        "origen": "orden de $a$",
        "seguridad": "Nula, y además las ternas caen en 15 planos.",
        "usos": "Ninguno. Solo interés histórico y didáctico.",
    },
    {
        "nombre": "Cuadrados medios",
        "generar": lambda n: ms.uniformes(n, SEMILLAS["cuadrados_medios"], DIGITOS),
        "periodo": "cola $326\\,565$ + ciclo $2\\,500$",
        "origen": "\\textbf{medido}",
        "seguridad": "Nula. Degenera: cinco puntos fijos con 4 dígitos.",
        "usos": "Ninguno. Obsoleto desde los años 50.",
    },
    {
        "nombre": "Mersenne Twister",
        "generar": lambda n: mersenne.uniformes(n, SEMILLAS["mersenne"]),
        "periodo": "$2^{19937}-1$",
        "origen": "especificación",
        "seguridad": "Nula: con 624 salidas se reconstruye el estado.",
        "usos": "Estándar de simulación. Motor de \\texttt{random} de Python.",
    },
    {
        "nombre": "Blum Blum Shub",
        "generar": lambda n: bbs.uniformes(n, SEMILLAS["bbs"]),
        "periodo": None,   # se calcula abajo
        "origen": "$\\mathrm{ord}_{p'q'}(2)$",
        "seguridad": "\\textbf{Criptográfica}: equivale a la residuosidad "
                     "cuadrática módulo $M$.",
        "usos": "Criptografía, donde la lentitud es aceptable.",
    },
]


def main():
    periodo_bbs, _ = bbs.periodo_estados()
    GENERADORES[-1]["periodo"] = f"$\\approx 2^{{{periodo_bbs.bit_length() - 1}}}$"

    print(f"\n--- Tabla comparativa (velocidad sobre n = {N_VELOCIDAD:,}) ---\n")
    print(f"  {'Generador':<22} {'s / 10^6':>10}  {'Período':<30} {'Origen'}")
    print("  " + "-" * 82)

    filas = []
    referencia = None
    for g in GENERADORES:
        segundos = cronometrar(g["generar"])
        if referencia is None:
            referencia = segundos
        relativo = segundos / referencia
        limpio = g["periodo"].replace("$", "").replace("\\,", " ").replace("\\approx ", "~")
        print(f"  {g['nombre']:<22} {segundos:>10.2f}  {limpio:<30} {g['origen']}")
        filas.append([g["nombre"], f"{segundos:.2f}", f"{relativo:.1f}$\\times$",
                      g["periodo"], g["seguridad"], g["usos"]])

    ruta = salida.guardar_tabla(
        ["Algoritmo", "s / $10^6$", "Rel.", "Período", "Seguridad", "Casos de uso"],
        filas, "tabla_comparativa.tex", alineacion="lrrlp{4cm}p{4cm}")
    print(f"\n  -> {ruta.name}")
    print()


if __name__ == "__main__":
    main()
