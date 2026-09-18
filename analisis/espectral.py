"""
Prueba espectral: las ternas (u_i, u_i+1, u_i+2) dentro del cubo [0,1)^3.

Por el teorema de estructura de red, si h = (h0, h1, h2) cumple
h0 + h1*a + h2*a^2 = 0 (mod m), entonces h0*u_i + h1*u_i+1 + h2*u_i+2 es entero
y todas las ternas caen en esa familia de planos. Para RANDU h = (9, -6, 1).
"""

import numpy as np

# Coeficientes del plano de RANDU: 9*u_i - 6*u_{i+1} + u_{i+2} = entero.
H_RANDU = np.array([9.0, -6.0, 1.0])


def ternas(u):
    """Ternas solapadas (u_i, u_{i+1}, u_{i+2}) como tres arreglos."""
    u = np.asarray(u, dtype=float)
    return u[:-2], u[1:-1], u[2:]


def combinacion(u, h=H_RANDU):
    """k_i = h0*u_i + h1*u_{i+1} + h2*u_{i+2}. Entero exacto si h cumple la red."""
    x, y, z = ternas(u)
    return h[0] * x + h[1] * y + h[2] * z


def contar_planos(k, tol=1e-9):
    """
    Planos ocupados: valores enteros distintos que toma k.

    Devuelve (cantidad, enteros_ordenados, desviacion_maxima), donde la
    desviacion maxima respecto al entero mas cercano mide que tan exacta es
    la relacion. Si es del orden de 1e-16 la identidad es algebraica; si es
    grande, k no es entero y no hay estructura de planos.
    """
    cercanos = np.rint(k)
    desviacion = float(np.max(np.abs(k - cercanos)))
    valores = np.unique(cercanos[np.abs(k - cercanos) < tol]).astype(int)
    return len(valores), valores, desviacion


def angulo_de_canto(h=H_RANDU, paso=1.0):
    """
    Elevacion y azimut desde los que los planos se ven de canto.

    La camara de matplotlib apunta en la direccion
        d = (cos(e)cos(a), cos(e)sin(a), sin(e)).
    Los planos aparecen como lineas cuando se mira paralelo a ellos, es decir
    cuando d es perpendicular a la normal h: d . h = 0. Se busca en una rejilla
    el par (e, a) que minimiza |d . h|, de modo que el angulo queda DERIVADO y
    no ajustado a ojo.
    """
    # Se busca solo con elevacion positiva: d y -d dan el mismo plano de
    # vista, y mirar desde arriba es mas legible.
    mejor, mejor_par = np.inf, (30.0, -60.0)
    for e in np.arange(0, 90 + paso, paso):
        for a in np.arange(-180, 180 + paso, paso):
            er, ar = np.radians(e), np.radians(a)
            d = np.array([np.cos(er) * np.cos(ar),
                          np.cos(er) * np.sin(ar),
                          np.sin(er)])
            valor = abs(float(d @ h))
            if valor < mejor:
                mejor, mejor_par = valor, (float(e), float(a))
    return mejor_par, mejor
