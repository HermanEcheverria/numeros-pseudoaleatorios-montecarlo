"""
Reduccion de varianza: variables antiteticas y variables de control.

El error de Monte Carlo es sigma/sqrt(N). Subir N es caro; bajar sigma es
gratis si se explota alguna estructura del integrando. Las dos tecnicas de
aqui explotan estructuras distintas, y por eso fallan en casos distintos.
"""

import numpy as np

Z_95 = 1.959963984540054


def simple(f, a, b, u):
    """Monte Carlo crudo: media de f sobre los puntos transformados."""
    x = a + (b - a) * np.asarray(u, dtype=float)
    return (b - a) * f(x)


def antiteticas(f, a, b, u):
    """
    Variables antiteticas: a cada U se le empareja 1-U.

    Con N/2 uniformes se producen N evaluaciones, asi que el presupuesto de
    evaluaciones se mantiene. La varianza del promedio del par es

        Var[(Y + Y')/2] = (Var Y + Var Y' + 2 Cov(Y, Y')) / 4,

    de modo que la tecnica ayuda solo si Cov(Y, Y') < 0. Eso ocurre cuando f
    es monotona. Si f es simetrica respecto al centro del intervalo, entonces
    f(a+b-x) = f(x), la correlacion es +1 y NO se gana nada.

    Devuelve las N/2 medias de pares.
    """
    u = np.asarray(u, dtype=float)
    x = a + (b - a) * u
    x_espejo = a + b - x
    return (b - a) * (f(x) + f(x_espejo)) / 2.0


def control(f, a, b, u, g, media_g):
    """
    Variables de control: se resta una variable de media conocida.

        Y_c = Y - c (g(X) - E[g(X)])

    El coeficiente optimo es c* = Cov(Y, g)/Var(g), que anula la mayor parte
    de la varianza si g correlaciona con f. Con ese c*,

        Var(Y_c) = Var(Y) (1 - rho^2),

    asi que el factor de reduccion es 1/(1 - rho^2): la tecnica sirve en la
    medida en que g se parezca a f. Si rho = 0 no hace nada.

    `c` se estima de la misma muestra, lo que introduce un sesgo de orden
    1/N, despreciable frente al error O(N^{-1/2}).

    Devuelve (muestra_corregida, c_optimo, rho).
    """
    u = np.asarray(u, dtype=float)
    x = a + (b - a) * u
    y = (b - a) * f(x)
    gx = g(x)

    cov = np.cov(y, gx, ddof=1)
    c = cov[0, 1] / cov[1, 1]
    rho = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    return y - c * (gx - media_g), float(c), float(rho)


def resumen(muestra, z=Z_95):
    """Estimacion puntual, varianza muestral, error estandar e IC."""
    muestra = np.asarray(muestra, dtype=float)
    n = muestra.size
    media = float(muestra.mean())
    var = float(muestra.var(ddof=1))
    se = np.sqrt(var / n)
    return media, var, se, media - z * se, media + z * se
