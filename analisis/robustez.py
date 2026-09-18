"""
Estudio de robustez sobre multiples semillas (pruebas de dos niveles).

Primer nivel:  las 4 pruebas sobre una secuencia -> 4 p-valores.
Segundo nivel: se repite con K semillas y se analiza la DISTRIBUCION de esos
p-valores, con dos criterios tomados de NIST SP 800-22:

  (a) la proporcion de rechazos debe rondar alpha;
  (b) bajo H0 los p-valores son Unif(0,1), lo que se contrasta con chi-cuadrado
      sobre 10 bins.

Las K semillas NO pueden elegirse al azar ni consecutivas: en un LCG de periodo
completo todas viven en el mismo ciclo, y dos semillas cercanas dan secuencias
que se solapan casi por completo. Aqui se obtienen saltando exactamente n
posiciones con la forma cerrada, lo que las deja disjuntas por construccion.
"""

import numpy as np
from scipy import stats

from analisis import pruebas


def saltar(x0, pasos, a, c, m):
    """
    Estado tras `pasos` iteraciones, sin iterar.

    Compone f(x) = a*x + c consigo misma por exponenciacion binaria: si
    f_j(x) = A_j*x + C_j, entonces f_j(f_k(x)) = A_j*A_k*x + A_j*C_k + C_j.
    Costo O(log pasos) en lugar de O(pasos).
    """
    A, C = 1, 0
    Ab, Cb = a % m, c % m
    while pasos:
        if pasos & 1:
            A, C = (Ab * A) % m, (Ab * C + Cb) % m
        Ab, Cb = (Ab * Ab) % m, (Ab * Cb + Cb) % m
        pasos >>= 1
    return (A * x0 + C) % m


def semillas_disjuntas(x0, K, n, a, c, m):
    """K semillas espaciadas n posiciones: las secuencias no comparten valores."""
    if K * n > m:
        raise ValueError(f"K*n = {K * n} excede el periodo maximo m = {m}.")
    return [saltar(x0, j * n, a, c, m) for j in range(K)]


def estudio(generador, semillas, n, k_bins=50):
    """
    Corre la bateria sobre cada semilla.

    `generador` es un callable gen(n, semilla) -> lista de uniformes.
    Devuelve {nombre_prueba: array de K p-valores}.
    """
    acumulado = {nombre: [] for nombre in pruebas.PRUEBAS}
    for semilla in semillas:
        u = generador(n, semilla)
        for nombre, (_, p) in pruebas.bateria(u, k=k_bins).items():
            acumulado[nombre].append(p)
    return {nombre: np.array(ps) for nombre, ps in acumulado.items()}


def tasa_rechazo(pvalores, alpha=pruebas.ALPHA):
    """
    Proporcion de rechazos y banda de aceptacion de NIST SP 800-22.

    Bajo H0 la proporcion esperada es alpha, con desviacion
    sqrt(alpha*(1-alpha)/K). NIST acepta el intervalo de +-3 desviaciones.
    Devuelve (proporcion, limite_inferior, limite_superior, dentro).
    """
    K = len(pvalores)
    proporcion = float(np.mean(pvalores < alpha))
    sigma = np.sqrt(alpha * (1 - alpha) / K)
    lo, hi = max(0.0, alpha - 3 * sigma), alpha + 3 * sigma
    return proporcion, lo, hi, lo <= proporcion <= hi


def uniformidad_pvalores(pvalores, bins=10):
    """
    Segundo nivel: chi-cuadrado de los p-valores contra Unif(0,1).

    Bajo H0 los p-valores de una prueba continua son exactamente Unif(0,1).
    Que se amontonen delata un sesgo que una sola corrida no detecta.
    """
    observados, _ = np.histogram(pvalores, bins=bins, range=(0.0, 1.0))
    esperados = len(pvalores) / bins
    estadistico = float(np.sum((observados - esperados) ** 2 / esperados))
    return estadistico, float(stats.chi2.sf(estadistico, df=bins - 1))
