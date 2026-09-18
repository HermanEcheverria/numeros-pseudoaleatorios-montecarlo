"""
Metodo de los cuadrados medios (von Neumann, 1949).

Se eleva el estado al cuadrado, se rellena con ceros a la izquierda hasta
2*digitos posiciones y se extraen los `digitos` centrales.

A diferencia del LCG, NO tiene teoria de periodo: no existe un Hull-Dobell que
lo prediga. El periodo depende de la semilla y solo se conoce midiendolo, que
es lo que hace `periodo_y_cola`.
"""


def siguiente(x, digitos):
    """Un paso: cuadrado, relleno a 2*digitos, extraccion de los centrales."""
    cuadrado = str(x * x).zfill(2 * digitos)
    inicio = digitos // 2
    return int(cuadrado[inicio:inicio + digitos])


def estados(n, semilla, digitos=4):
    """Los n primeros estados enteros. Se descarta la semilla, como en el LCG."""
    if digitos % 2 != 0:
        raise ValueError(f"`digitos` debe ser par, se recibió {digitos}.")
    if not 0 <= semilla < 10**digitos:
        raise ValueError(f"La semilla debe tener a lo sumo {digitos} dígitos.")

    salida = []
    x = semilla
    for _ in range(n):
        x = siguiente(x, digitos)
        salida.append(x)
    return salida


def uniformes(n, semilla, digitos=4):
    """n valores en [0, 1), normalizando entre 10**digitos."""
    escala = 10**digitos
    return [x / escala for x in estados(n, semilla, digitos)]


def periodo_y_cola(semilla, digitos=4, limite=1_000_000):
    """
    Periodo y cola de la orbita, midiendolos.

    Devuelve (periodo, cola). La cola es el numero de estados iniciales que no
    se vuelven a visitar: existe porque la transicion no es inyectiva, a
    diferencia del LCG. Un periodo de 1 con estado 0 es el "mecanismo cero".
    """
    visto = {}
    x = semilla
    for i in range(limite):
        if x in visto:
            return i - visto[x], visto[x]
        visto[x] = i
        x = siguiente(x, digitos)
    return None, None


def censo(digitos=4):
    """
    Recorre TODAS las semillas de `digitos` digitos y mide su orbita.

    Devuelve {semilla: (periodo, cola)}. Con 4 digitos son 10,000 semillas,
    barrido exhaustivo: el periodo no se estima, se conoce.
    """
    return {s: periodo_y_cola(s, digitos) for s in range(10**digitos)}
