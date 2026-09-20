"""
Estimador de Monte Carlo para integrales y su intervalo de confianza.

    I = int_a^b f(x) dx = (b-a) E[f(U)],   U ~ Unif(a,b)
    I_N = (b-a)/N * sum f(U_i)

El estimador es insesgado y su error estandar decae como O(N^{-1/2}); la
derivacion esta en el documento. Aqui solo se implementa.
"""

import math

import numpy as np

Z_95 = 1.959963984540054   # cuantil 0.975 de la normal estandar


def integrar(f, a, b, u, z=Z_95):
    """
    Estima int_a^b f(x) dx con las uniformes `u` y devuelve su IC.

    `u` son valores en [0,1) producidos por alguno de los generadores del
    proyecto; aqui se transforman a Unif(a,b) con x = a + (b-a)u.

    Devuelve (I_N, error_estandar, ic_inferior, ic_superior).
    """
    u = np.asarray(u, dtype=float)
    N = u.size
    if N < 2:
        raise ValueError("Se necesitan al menos 2 uniformes para estimar S_f.")

    fx = f(a + (b - a) * u)
    I = (b - a) * fx.mean()

    # ddof=1 (correccion de Bessel): np.std divide entre N por defecto, y la
    # derivacion del documento usa el estimador insesgado de la varianza.
    se = (b - a) * fx.std(ddof=1) / math.sqrt(N)
    return I, se, I - z * se, I + z * se


def cobertura(f, a, b, u, bloques, verdadero, z=Z_95):
    """
    Parte la secuencia en bloques disjuntos y cuenta cuantos IC atrapan el
    valor verdadero. Deberia rondar el 95%.

    Se usan bloques de una misma corrida y no semillas distintas porque en un
    LCG dos semillas dan secuencias relacionadas de forma determinista; los
    bloques disjuntos garantizan que ningun numero se reutilice.
    """
    u = np.asarray(u, dtype=float)
    tam = u.size // bloques
    intervalos = [integrar(f, a, b, u[k * tam:(k + 1) * tam], z)
                  for k in range(bloques)]
    aciertos = sum(lo <= verdadero <= hi for _, _, lo, hi in intervalos)
    return aciertos, intervalos
