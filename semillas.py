"""Registro central de semillas"""

N = 10_000

# Semillas de produccion: alimentan histogramas, pruebas y figuras.
# La del LCG es arbitraria a proposito: Hull-Dobell garantiza periodo completo
# para toda semilla, asi que el valor no cambia la calidad de la secuencia.
# La de RANDU no puede ser arbitraria: debe ser impar para alcanzar 2**29.
SEMILLAS = {
    "lcg": 25062003, # Fecha de cumpleaños mia 
    "randu": 111111111, # impar y suficientemente grande
}

# Semillas que ROMPEN el generador. Solo se usan en el laboratorio, para
# exhibir cada falla; nunca para producir resultados del informe.
SEMILLAS_DEGENERADAS = {
    "lcg_cero": 0,
    "randu_par": 2,
    "randu_constante": 2**30,
}
