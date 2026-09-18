# Números Pseudoaleatorios y Métodos de Monte Carlo

**Tarea 2 — Curso Libre de Configuración 2: Análisis de Datos**
Facultad de Ingeniería, Universidad del Istmo de Guatemala

**Autor:** Herman Andrés Echeverría Rojas — `hecheverria@unis.edu.gt`

---

## Requisitos

- Python 3.10 o superior

```bash
py -m pip install -r requirements.txt
```

## Cómo ejecutar

Cada script de `scripts/` reproduce una sección del informe: genera sus figuras
en `figuras/` y sus tablas en `resultados/`, siempre a partir de las semillas
registradas en `semillas.py`.

```bash
py -m scripts.<nombre_del_script>
```

## Estructura

| Carpeta | Contenido |
|---|---|
| `prngs/` | Los cinco generadores, implementados desde cero |
| `analisis/` | Pruebas de hipótesis, prueba espectral y utilidades de figuras |
| `montecarlo/` | Estimadores de Monte Carlo y técnicas de reducción de varianza |
| `scripts/` | Un script por sección; produce las figuras y tablas del informe |
| `figuras/` | Figuras generadas (no se editan a mano) |
| `resultados/` | Tablas `.tex` generadas (no se editan a mano) |

## Reproducibilidad

Todas las semillas viven en `semillas.py`, en un solo lugar. Cada generador
recibe su semilla como argumento y es completamente determinista: **misma
semilla ⇒ misma secuencia**.

Ningún número del informe está transcrito a mano: los scripts escriben las
tablas como archivos `.tex` en `resultados/`, listas para que el documento las
incorpore con `\input{}`.

## Restricción del enunciado

No se usa `numpy.random` ni `random` como motor de ningún generador. Todos los
algoritmos están implementados a partir de su recurrencia matemática.

`numpy` se usa solo para manejo de arreglos y aritmética vectorizada,
`scipy.stats` para obtener p-valores de distribuciones teóricas, y `matplotlib`
para graficar. `random` y `secrets` aparecen únicamente en la comparación que
pide el inciso (d) de la Parte 1.

## Estado de avance

- [x] Fundamento matemático común (periodicidad, forma cerrada, Hull–Dobell, estructura de red)
- [ ] LCG
- [ ] RANDU
- [ ] Cuadrados medios (von Neumann)
- [ ] Blum Blum Shub
- [ ] Mersenne Twister (MT19937)
- [ ] Prueba espectral visual
- [ ] Comparación con `random` y `secrets`
- [ ] Monte Carlo: fundamento del estimador y su error
- [ ] Monte Carlo: integrales en una dimensión y convergencia
- [ ] Monte Carlo: alta dimensión y cuadratura en rejilla
- [ ] Monte Carlo: reducción de varianza
- [ ] Monte Carlo: el caso RANDU

## Referencias

- Hull, T. E. y Dobell, A. R. (1962). *Random Number Generators*. SIAM Review, 4(3), 230–254.
- Marsaglia, G. (1968). *Random Numbers Fall Mainly in the Planes*. PNAS, 61(1), 25–28.
- Knuth, D. E. (1997). *The Art of Computer Programming, Vol. 2*, 3.ª ed., cap. 3.
- Matsumoto, M. y Nishimura, T. (1998). *Mersenne Twister*. ACM TOMACS, 8(1), 3–30.
- Blum, L., Blum, M. y Shub, M. (1986). *A Simple Unpredictable Pseudo-Random Number Generator*. SIAM J. Comput., 15(2), 364–383.
