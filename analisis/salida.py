"""
Todo lo que el proyecto escribe a disco para el informe: figuras y tablas.

Las tablas se guardan como .tex en resultados/ para que el documento las
incorpore con \\input{}: ningun numero se transcribe a mano.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAIZ = Path(__file__).resolve().parent.parent
FIGURAS = RAIZ / "figuras"
RESULTADOS = RAIZ / "resultados"

plt.rcParams.update({
    "figure.dpi": 150,
    "font.size": 9,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "axes.axisbelow": True,
})


def _guardar(fig, nombre):
    FIGURAS.mkdir(exist_ok=True)
    ruta = FIGURAS / nombre
    fig.savefig(ruta, bbox_inches="tight")
    plt.close(fig)
    return ruta


def histogramas(secuencias, nombre_archivo, bins=50, titulo_extra=""):
    """
    Histogramas normalizados lado a lado, con la densidad teorica en y = 1.

    `secuencias` es una lista de (etiqueta, valores, semilla).
    """
    fig, ejes = plt.subplots(1, len(secuencias),
                             figsize=(4.2 * len(secuencias), 3.2), sharey=True)
    if len(secuencias) == 1:
        ejes = [ejes]

    for eje, (etiqueta, u, semilla) in zip(ejes, secuencias):
        eje.hist(u, bins=bins, range=(0, 1), density=True,
                 color="#4C72B0", edgecolor="white", linewidth=0.4)
        eje.axhline(1.0, color="#C44E52", linestyle="--", linewidth=1.2,
                    label="densidad Unif(0,1)")
        eje.set_title(f"{etiqueta}\nsemilla = {semilla}, n = {len(u):,}")
        eje.set_xlabel("$u$")
        eje.set_xlim(0, 1)
        eje.legend(fontsize=7, loc="lower right")

    ejes[0].set_ylabel("densidad")
    if titulo_extra:
        fig.suptitle(titulo_extra)
    return _guardar(fig, nombre_archivo)


def histogramas_pvalores(estudios, nombre_archivo, bins=10):
    """
    Rejilla de histogramas de p-valores: una fila por generador, una columna
    por prueba. Bajo H0 cada uno deberia ser plano en 1.
    """
    generadores = list(estudios)
    nombres_pruebas = list(estudios[generadores[0]])

    fig, ejes = plt.subplots(len(generadores), len(nombres_pruebas),
                             figsize=(2.6 * len(nombres_pruebas),
                                      2.3 * len(generadores)),
                             squeeze=False, sharex=True)

    for i, generador in enumerate(generadores):
        for j, prueba in enumerate(nombres_pruebas):
            eje = ejes[i][j]
            pvals = estudios[generador][prueba]
            eje.hist(pvals, bins=bins, range=(0, 1), density=True,
                     color="#4C72B0", edgecolor="white", linewidth=0.4)
            eje.axhline(1.0, color="#C44E52", linestyle="--", linewidth=1.0)
            eje.set_xlim(0, 1)
            if i == 0:
                eje.set_title(prueba, fontsize=9)
            if j == 0:
                eje.set_ylabel(f"{generador}\ndensidad", fontsize=8)
            if i == len(generadores) - 1:
                eje.set_xlabel("$p$-valor", fontsize=8)

    return _guardar(fig, nombre_archivo)


def guardar_tabla(encabezados, filas, nombre_archivo, alineacion=None):
    """Escribe un tabular de LaTeX en resultados/, listo para \\input{}."""
    RESULTADOS.mkdir(exist_ok=True)
    if alineacion is None:
        alineacion = "l" + "r" * (len(encabezados) - 1)

    lineas = [f"\\begin{{tabular}}{{{alineacion}}}", "\\toprule",
              " & ".join(encabezados) + " \\\\", "\\midrule"]
    lineas += [" & ".join(str(c) for c in fila) + " \\\\" for fila in filas]
    lineas += ["\\bottomrule", "\\end{tabular}"]

    ruta = RESULTADOS / nombre_archivo
    ruta.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    return ruta


def formato_p(p):
    """p-valores: notacion cientifica cuando son diminutos, para no ver 0.0000."""
    return f"{p:.4f}" if p >= 1e-4 else f"{p:.2e}"
