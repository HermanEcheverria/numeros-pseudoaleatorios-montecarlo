"""Generador congruencial lineal y RANDU: x_{n+1} = (a*x_n + c) mod m."""

from math import gcd

PARAMS_NUMERICAL_RECIPES = {"a": 1664525, "c": 1013904223, "m": 2**32}
PARAMS_GLIBC = {"a": 1103515245, "c": 12345, "m": 2**31}
PARAMS_RANDU = {"a": 65539, "c": 0, "m": 2**31}


def estados(n, semilla, a, c, m, validar=True):
    """
    Los n primeros estados enteros x_1, ..., x_n.

    Se descarta x_0: la primera salida es x_1, no la semilla, que es un valor
    elegido por el usuario y no producido por el generador.

    `validar=False` permite exhibir los casos degenerados en el laboratorio.
    """
    if n <= 0:
        raise ValueError(f"n debe ser positivo, se recibió {n}.")
    if m <= 1:
        raise ValueError(f"El módulo m debe ser mayor que 1, se recibió {m}.")
    if not 0 <= semilla < m:
        raise ValueError(f"La semilla debe estar en [0, {m}), se recibió {semilla}.")
    if not 0 <= a < m or not 0 <= c < m:
        raise ValueError(f"'a' y 'c' deben estar en [0, {m}).")

    if validar:
        if a == 0:
            raise ValueError("a = 0: la secuencia es la constante c/m.")
        if a == 1:
            raise ValueError("a = 1: la secuencia es una progresión aritmética.")
        if c == 0 and semilla == 0:
            raise ValueError("c = 0 con semilla 0: el cero es punto fijo.")

    salida = []
    x = semilla
    for _ in range(n):
        x = (a * x + c) % m
        salida.append(x)
    return salida


def uniformes(n, semilla, a, c, m, validar=True):
    """
    n uniformes u_n = x_n / m en [0, 1).

    Se divide entre m y no entre m-1 para que el máximo sea (m-1)/m < 1: con
    m-1 la salida u = 1 rompería cualquier transformación posterior que use
    log(1-u) o la inversa de una acumulada.
    """
    return [x / m for x in estados(n, semilla, a, c, m, validar=validar)]


def randu(n, semilla, validar=True):
    """RANDU (a=65539, c=0, m=2**31). Necesita semilla impar para su período máximo."""
    return uniformes(n, semilla, validar=validar, **PARAMS_RANDU)


def cumple_hull_dobell(a, c, m):
    """
    Verifica las tres condiciones de Hull-Dobell para período completo.

    Devuelve (cumple_todas, detalle), donde detalle lista cada condición con
    su veredicto para poder reportar cuál falló.
    """
    cond1 = gcd(c, m) == 1
    factores = _factores_primos(m)
    cond2 = all((a - 1) % q == 0 for q in factores)
    cond3 = (a - 1) % 4 == 0 if m % 4 == 0 else True

    detalle = [
        (f"gcd(c, m) = gcd({c}, {m}) = {gcd(c, m)}, debe ser 1", cond1),
        (f"a - 1 = {a - 1} divisible por los primos de m {sorted(factores)}", cond2),
        (f"4 | m es {m % 4 == 0}; si lo es, 4 | (a-1) = {(a - 1) % 4 == 0}", cond3),
    ]
    return all([cond1, cond2, cond3]), detalle


def _factores_primos(n):
    """Factores primos distintos de n. Suficiente para m <= 2**32."""
    factores = set()
    d = 2
    while d * d <= n:
        while n % d == 0:
            factores.add(d)
            n //= d
        d += 1
    if n > 1:
        factores.add(n)
    return factores


def periodo(semilla, a, c, m, limite=None):
    """
    Período y cola de la órbita que arranca en `semilla`, por fuerza bruta.

    Devuelve (periodo, cola). Costo O(periodo): solo apto para m pequeño.
    Para m grande usar `orden_multiplicativo`.
    """
    if limite is None:
        limite = m + 1
    visto = {}
    x = semilla
    for i in range(limite):
        if x in visto:
            return i - visto[x], visto[x]
        visto[x] = i
        x = (a * x + c) % m
    return None, None


def orden_multiplicativo(a, m):
    """
    Menor t con a**t = 1 (mod m), para m potencia de 2.

    En un LCG multiplicativo con semilla impar ese orden es el período. Se
    prueban potencias de 2 con pow() modular, que es O(log t): iterar la
    recurrencia 2**29 veces sería inviable.
    """
    if gcd(a, m) != 1:
        raise ValueError(f"a = {a} no es invertible módulo {m}.")
    t = 1
    while pow(a, t, m) != 1:
        t *= 2
        if t > m:
            raise ValueError("No se encontró el orden; m no es potencia de 2.")
    return t


def verificar_identidad_randu(n, semilla):
    """
    Comprueba x_{n+2} = 6*x_{n+1} - 9*x_n (mod 2**31) en todas las ternas.

    Devuelve (se_cumple_siempre, ternas_revisadas). Se trabaja sobre los
    estados enteros porque la identidad es exacta: con los valores ya
    normalizados el redondeo de punto flotante daría falsos negativos.
    """
    m = PARAMS_RANDU["m"]
    xs = estados(n, semilla, **PARAMS_RANDU)
    ternas = len(xs) - 2
    for i in range(ternas):
        if (6 * xs[i + 1] - 9 * xs[i]) % m != xs[i + 2]:
            return False, ternas
    return True, ternas
