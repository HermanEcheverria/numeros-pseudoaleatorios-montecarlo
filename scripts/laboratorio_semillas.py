"""
Laboratorio de semillas: dónde se rompe el LCG.

Comprueba empíricamente los resultados teóricos del fundamento. La
interpretación de cada hallazgo va en el documento, no aquí.

    py -m scripts.laboratorio_semillas
"""

from prngs import lcg
from semillas import SEMILLAS, SEMILLAS_DEGENERADAS


def titulo(texto):
    print(f"\n{'=' * 72}\n{texto}\n{'=' * 72}")


def hull_dobell():
    titulo("1. Hull-Dobell sobre los tres juegos de parámetros")

    for nombre, p in [("Numerical Recipes", lcg.PARAMS_NUMERICAL_RECIPES),
                      ("glibc", lcg.PARAMS_GLIBC),
                      ("RANDU", lcg.PARAMS_RANDU)]:
        cumple, detalle = lcg.cumple_hull_dobell(p["a"], p["c"], p["m"])
        print(f"\n{nombre}:  a = {p['a']}, c = {p['c']}, "
              f"m = 2**{p['m'].bit_length() - 1}")
        for descripcion, ok in detalle:
            print(f"    [{'OK ' if ok else 'NO'}] {descripcion}")
        print("    =>", f"período completo = {p['m']}" if cumple
              else "NO alcanza período completo")


def verificacion_modulo_pequeno():
    titulo("2. Períodos por fuerza bruta sobre las 256 semillas de m = 2**8")

    m = 2**8
    print(f"\n  {'a':>4} {'c':>4}  {'HD':>3}   períodos observados")
    print("  " + "-" * 60)
    for a, c in [(13, 7), (13, 8), (11, 7)]:
        cumple, _ = lcg.cumple_hull_dobell(a, c, m)
        periodos = sorted({lcg.periodo(s, a, c, m)[0] for s in range(m)})
        print(f"  {a:>4} {c:>4}  {'sí' if cumple else 'no':>3}   {periodos}")


def casos_degenerados():
    titulo("3. Casos degenerados")

    m = 2**31
    for etiqueta, params, semilla in [
        ("c = 0, semilla = 0", {"a": 65539, "c": 0},
         SEMILLAS_DEGENERADAS["lcg_cero"]),
        ("a = 0",              {"a": 0, "c": 7},     12345),
        ("a = 1",              {"a": 1, "c": 1},     0),
    ]:
        u = lcg.uniformes(6, semilla, m=m, validar=False, **params)
        print(f"\n  {etiqueta:<20} {[f'{v:.4e}' for v in u]}")


def randu_semillas():
    titulo("4. RANDU: el período contra los factores 2 de la semilla")

    a, m = lcg.PARAMS_RANDU["a"], lcg.PARAMS_RANDU["m"]
    orden = lcg.orden_multiplicativo(a, m)
    print(f"\n  Semilla impar -> período = ord_(2**31)(65539) = {orden} "
          f"= 2**{orden.bit_length() - 1}")

    print(f"\n  {'semilla':>16} {'= 2**s':>10} {'período':>14}  {'distintas de 64':>17}")
    print("  " + "-" * 62)
    for s in [0, 1, 2, 8, 16, 26, 27, 28, 29, 30]:
        # Un factor 2**s en la semilla hace que el generador se comporte como
        # si el módulo fuera 2**(31-s).
        m_efectivo = 2 ** (31 - s)
        p = 1 if m_efectivo <= 2 else lcg.orden_multiplicativo(a % m_efectivo,
                                                              m_efectivo)
        distintas = len(set(lcg.uniformes(64, 2**s, a=a, c=0, m=m, validar=False)))
        print(f"  {2**s:>16} {'2**' + str(s):>10} {p:>14}  {distintas:>17}")

    for clave in ("randu_par", "randu_constante"):
        semilla = SEMILLAS_DEGENERADAS[clave]
        u = lcg.uniformes(6, semilla, a=a, c=0, m=m, validar=False)
        print(f"\n  {clave} = {semilla}: {[f'{v:.4f}' for v in u]}")
    print(f"  Semilla del informe: {SEMILLAS['randu']} "
          f"(impar: {SEMILLAS['randu'] % 2 == 1})")


def identidad_randu():
    titulo("5. Identidad x_(n+2) = 6*x_(n+1) - 9*x_n (mod 2**31)")

    ok, ternas = lcg.verificar_identidad_randu(10_000, SEMILLAS["randu"])
    print(f"\n  Semilla {SEMILLAS['randu']}: se cumple en las "
          f"{ternas:,} ternas = {ok}")


def main():
    hull_dobell()
    verificacion_modulo_pequeno()
    casos_degenerados()
    randu_semillas()
    identidad_randu()
    print()


if __name__ == "__main__":
    main()
