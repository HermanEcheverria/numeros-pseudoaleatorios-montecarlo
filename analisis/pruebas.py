"""
Bateria de pruebas de hipotesis sobre una secuencia de uniformes.

Uniformidad:   chi-cuadrado y Kolmogorov-Smirnov
Independencia: autocorrelacion de retardo 1 y rachas

Con alpha = 0.05 se rechaza H0 si p < alpha. Aqui se BUSCA no rechazar, pero
"no se rechaza H0" no equivale a "se demuestra H0".
"""

import numpy as np
from scipy import stats

ALPHA = 0.05


def chi2_uniformidad(u, k=50):
    """
    Bondad de ajuste contra Unif(0,1) con k bins iguales.

        chi2 = sum (O_i - E_i)^2 / E_i,   E_i = n/k,   gl = k - 1

    Se pierde un grado de libertad porque los O_i estan obligados a sumar n.
    """
    u = np.asarray(u)
    n = len(u)
    esperados = n / k
    if esperados < 5:
        raise ValueError(f"Frecuencia esperada {esperados:.1f} < 5; reducir k.")

    # range=(0,1) fija los limites: sin el, numpy usaria el min/max de los
    # datos y no se estaria midiendo uniformidad sobre [0,1).
    observados, _ = np.histogram(u, bins=k, range=(0.0, 1.0))
    if observados.sum() != n:
        raise ValueError(f"{n - observados.sum()} valores fuera de [0, 1).")

    estadistico = np.sum((observados - esperados) ** 2 / esperados)
    return estadistico, stats.chi2.sf(estadistico, df=k - 1), k - 1


def ks_uniformidad(u):
    """
    Kolmogorov-Smirnov de una muestra contra Unif(0,1).

        D = sup_x |F_n(x) - x|

    No requiere elegir bins, a diferencia de chi-cuadrado.
    """
    r = stats.kstest(np.asarray(u), "uniform")
    return r.statistic, r.pvalue


def autocorrelacion_lag1(u):
    """
    Correlacion de la serie contra si misma desplazada un lugar.

    Bajo H0 se espera r1 ~ 0, y z = r1*sqrt(n) es aproximadamente N(0,1).
    El p-valor es a dos colas: tanto la tendencia como la alternancia son
    evidencia de dependencia.
    """
    u = np.asarray(u)
    n = len(u)
    media = u.mean()
    r1 = np.sum((u[:-1] - media) * (u[1:] - media)) / np.sum((u - media) ** 2)
    z = r1 * np.sqrt(n)
    return r1, z, 2 * stats.norm.sf(abs(z))


def rachas(u):
    """
    Prueba de rachas respecto a la mediana.

    Se marca + si u_i >= mediana y - si no; R es el numero de bloques
    consecutivos del mismo simbolo. La variante de la mediana garantiza
    n1 ~ n2 por construccion, lo que evita que un desbalance distorsione
    E[R] y Var[R].
    """
    u = np.asarray(u)
    n = len(u)
    binario = (u >= np.median(u)).astype(int)
    R = int(np.sum(np.diff(binario) != 0) + 1)

    # float() evita el desbordamiento silencioso de int64: el numerador de
    # Var_R crece como n^4/4 y para n >= ~75,000 excede el maximo de int64.
    n1 = float(np.sum(binario == 1))
    n2 = float(np.sum(binario == 0))

    E_R = 2 * n1 * n2 / n + 1
    Var_R = 2 * n1 * n2 * (2 * n1 * n2 - n) / (n**2 * (n - 1))
    z = (R - E_R) / np.sqrt(Var_R)
    return R, z, 2 * stats.norm.sf(abs(z))


PRUEBAS = ["chi2", "ks", "autocorrelacion", "rachas"]


def bateria(u, k=50):
    """Las cuatro pruebas sobre u. Devuelve {nombre: (estadistico, p_valor)}."""
    chi_e, chi_p, _ = chi2_uniformidad(u, k=k)
    ks_d, ks_p = ks_uniformidad(u)
    r1, _, auto_p = autocorrelacion_lag1(u)
    R, _, rachas_p = rachas(u)
    return {
        "chi2": (chi_e, chi_p),
        "ks": (ks_d, ks_p),
        "autocorrelacion": (r1, auto_p),
        "rachas": (float(R), rachas_p),
    }


def decision(p, alpha=ALPHA):
    return "no se rechaza" if p >= alpha else "SE RECHAZA"
