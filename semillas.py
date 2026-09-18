"""Registro central de semillas"""

N = 10_000

# Semillas de produccion: alimentan histogramas, pruebas y figuras.
# La del LCG es arbitraria a proposito: Hull-Dobell garantiza periodo completo
# para toda semilla, asi que el valor no cambia la calidad de la secuencia.
# La de RANDU no puede ser arbitraria: debe ser impar para alcanzar 2**29.
SEMILLAS = {
    "lcg": 25062003, # Fecha de cumpleaños mia
    "randu": 111111111, # impar y suficientemente grande
    "cuadrados_medios": 250620032506, # 12 digitos: cola 326,565 y periodo 2,500
    "mersenne": 25062003, # cualquier valor de 32 bits sirve
}

# Semilla de referencia de la especificacion MT19937, para verificar la
# implementacion contra el vector publicado (OEIS A221557).
SEMILLA_REFERENCIA_MT = 5489

# Cuadrados medios no tiene teoria de periodo: el numero de digitos decide
# cuanto dura la orbita y hay que medirlo. Con 4 digitos el periodo maximo
# sobre las 10,000 semillas posibles es 4, inservible para n = 10,000.
DIGITOS_CUADRADOS_MEDIOS = 12

# Semillas que ROMPEN el generador. Solo se usan en el laboratorio, para
# exhibir cada falla; nunca para producir resultados del informe.
SEMILLAS_DEGENERADAS = {
    "lcg_cero": 0,
    "randu_par": 2,
    "randu_constante": 2**30,
}
