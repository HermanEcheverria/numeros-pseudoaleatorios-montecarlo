"""
Volumen de la bola unitaria en d dimensiones: Monte Carlo contra cuadratura.

    V_d = Vol{ x en R^d : ||x||_2 <= 1 } = pi^(d/2) / Gamma(d/2 + 1)

El valor cerrado se usa SOLO como verdad de referencia, nunca dentro de los
estimadores.
"""

import math

import numpy as np

Z_95 = 1.959963984540054


def volumen_exacto(d):
    """pi^(d/2) / Gamma(d/2 + 1)."""
    return math.pi ** (d / 2) / math.gamma(d / 2 + 1)


def monte_carlo(d, u, z=Z_95):
    """
    Estima V_d muestreando en el cubo [-1,1]^d y contando los puntos dentro.

    `u` son N*d uniformes en [0,1), que se reacomodan en N puntos de d
    coordenadas y se llevan a [-1,1] con x = 2u - 1.

    La fraccion p que cae dentro es una proporcion binomial, asi que
    Var(p) = p(1-p)/N y el IC sale de ahi, multiplicado por el volumen del
    cubo 2^d.

    ATENCION: si ningun punto cae dentro, p = 0 y el IC normal colapsa a
    [0, 0], que no contiene el valor real. No es un error de calculo: es que
    la aproximacion normal exige varios aciertos. En ese caso se usa la regla
    de tres, que da como cota superior al 95% el valor 3/N.

    Devuelve (V_estimado, error_estandar, ic_inferior, ic_superior, N, aciertos).
    """
    u = np.asarray(u, dtype=float)
    N = u.size // d
    puntos = 2.0 * u[:N * d].reshape(N, d) - 1.0
    dentro = (np.sum(puntos * puntos, axis=1) <= 1.0)

    aciertos = int(dentro.sum())
    p = aciertos / N
    volumen_cubo = 2.0**d
    V = volumen_cubo * p

    if aciertos == 0:
        # Regla de tres: con 0 exitos en N ensayos, el limite superior del
        # IC del 95% para la proporcion es aproximadamente 3/N.
        return V, float("nan"), 0.0, volumen_cubo * 3.0 / N, N, 0

    se = volumen_cubo * math.sqrt(p * (1 - p) / N)
    return V, se, V - z * se, V + z * se, N, aciertos


def cuadratura_rejilla(d, m):
    """
    Estimacion por rejilla: m subintervalos por eje, integrando evaluado en el
    centro de cada celda. Requiere m^d evaluaciones.

    Solo es viable para d pequeno; para d grande se devuelve None porque m^d
    excede cualquier presupuesto. Esa imposibilidad ES el resultado.
    """
    if m**d > 5_000_000:
        return None, m**d

    centros = (np.arange(m) + 0.5) / m * 2.0 - 1.0      # centros en [-1,1]
    rejilla = np.meshgrid(*([centros] * d), indexing="ij")
    cuadrados = sum(coord**2 for coord in rejilla)
    fraccion = float(np.mean(cuadrados <= 1.0))
    return 2.0**d * fraccion, m**d


def puntos_para_rejilla(d, m=10):
    """Evaluaciones que necesitaria una rejilla de m subdivisiones por eje."""
    return m**d
