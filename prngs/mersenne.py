"""
Mersenne Twister MT19937 (Matsumoto y Nishimura, 1998).

Recurrencia lineal sobre GF(2) con 624 palabras de 32 bits de estado. Tres
etapas: inicializacion desde la semilla, torcedura (twist) del estado completo
y templado (tempering) de cada palabra antes de emitirla.

Periodo 2**19937 - 1, primo de Mersenne, de donde viene el nombre.

Las constantes son las de la especificacion MT19937; la implementacion se
verifica contra el vector de referencia de la semilla 5489 (OEIS A221557).
"""

# Parametros de MT19937
N_ESTADO = 624
M_DESPLAZAMIENTO = 397
A_TORCEDURA = 0x9908B0DF
F_INICIAL = 1812433253

MASCARA_32 = 0xFFFFFFFF
MASCARA_SUPERIOR = 0x80000000   # bit mas significativo
MASCARA_INFERIOR = 0x7FFFFFFF   # los 31 bits restantes

# Constantes del templado
U, S, T, L = 11, 7, 15, 18
B_TEMPLADO = 0x9D2C5680
C_TEMPLADO = 0xEFC60000


class MersenneTwister:
    """Estado de 624 palabras. Una instancia por semilla."""

    def __init__(self, semilla):
        if not 0 <= semilla <= MASCARA_32:
            raise ValueError(f"La semilla debe caber en 32 bits, se recibió {semilla}.")
        self.mt = [0] * N_ESTADO
        self.mt[0] = semilla
        for i in range(1, N_ESTADO):
            previo = self.mt[i - 1] ^ (self.mt[i - 1] >> 30)
            self.mt[i] = (F_INICIAL * previo + i) & MASCARA_32
        self.indice = N_ESTADO   # fuerza una torcedura en la primera extraccion

    def _torcer(self):
        """Regenera las 624 palabras de golpe."""
        for i in range(N_ESTADO):
            # Se concatena el bit alto de mt[i] con los 31 bits bajos del siguiente
            x = ((self.mt[i] & MASCARA_SUPERIOR) |
                 (self.mt[(i + 1) % N_ESTADO] & MASCARA_INFERIOR))
            xa = x >> 1
            if x % 2:                       # multiplicacion por la matriz A
                xa ^= A_TORCEDURA
            self.mt[i] = self.mt[(i + M_DESPLAZAMIENTO) % N_ESTADO] ^ xa
        self.indice = 0

    def extraer(self):
        """Siguiente palabra de 32 bits, ya templada."""
        if self.indice >= N_ESTADO:
            self._torcer()

        y = self.mt[self.indice]
        self.indice += 1

        # Templado: mejora la equidistribucion de los bits. Es invertible, asi
        # que no aporta seguridad; con 624 salidas se recupera el estado.
        y ^= y >> U
        y ^= (y << S) & B_TEMPLADO
        y ^= (y << T) & C_TEMPLADO
        y ^= y >> L
        return y & MASCARA_32


def estados(n, semilla):
    """Las n primeras palabras de 32 bits."""
    generador = MersenneTwister(semilla)
    return [generador.extraer() for _ in range(n)]


def uniformes(n, semilla):
    """n valores en [0, 1), dividiendo entre 2**32."""
    return [x / 2**32 for x in estados(n, semilla)]
