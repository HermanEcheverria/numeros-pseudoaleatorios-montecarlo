"""
Blum Blum Shub (Blum, Blum y Shub, 1986).

    x_{n+1} = x_n^2 mod M,   M = p*q

con p y q primos de Blum (congruentes con 3 modulo 4). Es el unico de los cinco
generadores criptograficamente seguro: distinguir su salida de una secuencia
aleatoria es tan dificil como resolver el problema de la residuosidad
cuadratica modulo M, que a su vez se apoya en la dificultad de factorizar M.

Solo se emiten los log2(log2(M)) bits menos significativos de cada estado: es
la cantidad maxima para la que se conserva la garantia de seguridad. Emitir el
estado completo la destruiria, porque revelaria la orbita.
"""

from math import gcd, isqrt

# p y q son primos SEGUROS (p = 2p'+1 con p' primo) y de BLUM (p = 3 mod 4).
# Ambas propiedades se comprueban en `verificar_primos`.
P = 2147483783
Q = 2148484187
M = P * Q


def _es_primo(n):
    """Miller-Rabin determinista para n < 3.3e24 con estas bases."""
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def verificar_primos(p=P, q=Q):
    """
    Comprueba las condiciones que p y q deben cumplir.

    Devuelve una lista de (descripcion, se_cumple).
    """
    return [
        (f"p = {p} es primo", _es_primo(p)),
        (f"q = {q} es primo", _es_primo(q)),
        (f"p = 3 (mod 4): p mod 4 = {p % 4}", p % 4 == 3),
        (f"q = 3 (mod 4): q mod 4 = {q % 4}", q % 4 == 3),
        (f"p es seguro: (p-1)/2 = {(p - 1) // 2} primo", _es_primo((p - 1) // 2)),
        (f"q es seguro: (q-1)/2 = {(q - 1) // 2} primo", _es_primo((q - 1) // 2)),
        (f"gcd((p-3)/2, (q-3)/2) = {gcd((p - 3) // 2, (q - 3) // 2)} es pequeño",
         gcd((p - 3) // 2, (q - 3) // 2) <= 4),
    ]


def bits_por_paso(m=M):
    """log2(log2(M)) redondeado hacia abajo: los bits que es seguro emitir."""
    return m.bit_length().bit_length() - 1


def estado_inicial(semilla, m=M):
    """
    x_0 = semilla^2 mod M, como pide la construccion.

    Elevar al cuadrado garantiza que x_0 sea un residuo cuadratico, que es lo
    que hace que la orbita viva en el subgrupo correcto. La semilla debe ser
    coprima con M: si compartiera un factor, la orbita colapsaria en un
    subgrupo diminuto y ademas revelaria la factorizacion.
    """
    if gcd(semilla, m) != 1:
        raise ValueError(f"La semilla {semilla} no es coprima con M.")
    if semilla % m in (0, 1):
        raise ValueError("Semilla trivial: 0 y 1 son puntos fijos.")
    return semilla * semilla % m


def estados(n, semilla, p=P, q=Q):
    """Los n primeros estados enteros de la recurrencia."""
    m = p * q
    x = estado_inicial(semilla, m)
    salida = []
    for _ in range(n):
        x = x * x % m
        salida.append(x)
    return salida


def uniformes(n, semilla, p=P, q=Q, palabra=32):
    """
    n valores en [0, 1), armados con los bits bajos de estados sucesivos.

    Cada salida se compone de `palabra` bits, tomados de `bits_por_paso()` en
    `bits_por_paso()`. Por eso BBS necesita varias exponenciaciones modulares
    por numero, mientras que un LCG necesita una multiplicacion: ahi esta la
    diferencia de velocidad.
    """
    m = p * q
    k = bits_por_paso(m)
    mascara = (1 << k) - 1

    x = estado_inicial(semilla, m)
    salida, acumulado, disponibles = [], 0, 0
    while len(salida) < n:
        x = x * x % m
        acumulado = (acumulado << k) | (x & mascara)
        disponibles += k
        if disponibles >= palabra:
            sobran = disponibles - palabra
            salida.append((acumulado >> sobran) / 2**palabra)
            acumulado &= (1 << sobran) - 1
            disponibles = sobran
    return salida


def _factores_primos(n):
    """Factores primos distintos de n, por division tentativa."""
    factores = set()
    d = 2
    while d * d <= n:
        while n % d == 0:
            factores.add(d)
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        factores.add(n)
    return factores


def _orden_multiplicativo(a, n, exponente):
    """
    Menor t > 0 con a^t = 1 (mod n), sabiendo que t divide `exponente`.

    Se parte del exponente y se va dividiendo entre cada factor primo mientras
    la congruencia siga cumpliendose. Evita recorrer los divisores uno por uno.
    """
    if gcd(a, n) != 1:
        raise ValueError(f"{a} no es invertible módulo {n}.")
    orden = exponente
    for f in _factores_primos(exponente):
        while orden % f == 0 and pow(a, orden // f, n) == 1:
            orden //= f
    return orden


def periodo_estados(p=P, q=Q):
    """
    Periodo de la sucesion de estados, sin iterarla.

    Como x_n = x_0^(2^n) mod M, la sucesion se repite cuando 2^n vuelve a su
    valor inicial modulo el orden de x_0, asi que el periodo es el orden
    multiplicativo de 2 modulo ord(x_0).

    Con p y q seguros, lambda(M) = 2*p'*q' con p' = (p-1)/2 y q' = (q-1)/2
    primos. El orden NO se calcula modulo lambda(M): al ser par, 2 no es
    invertible ahi. Como x_0 es residuo cuadratico, su orden divide
    lambda(M)/2 = p'*q', que es impar, y el periodo es ord_(p'q')(2).

    Ese orden divide el exponente del grupo de unidades modulo p'*q', que es
    lcm(p'-1, q'-1), no p'*q'.
    """
    pp, qq = (p - 1) // 2, (q - 1) // 2

    # La formula del exponente vale solo si p' y q' son primos, es decir si p
    # y q son primos SEGUROS. Sin este guardian, con p' compuesto el exponente
    # sale mal y la funcion devuelve un periodo incorrecto sin avisar.
    if not (_es_primo(pp) and _es_primo(qq)):
        raise ValueError(
            f"p y q deben ser primos seguros: (p-1)/2 = {pp} y (q-1)/2 = {qq} "
            "tienen que ser primos."
        )

    lam = 2 * pp * qq
    n = pp * qq                                  # lambda(M)/2, impar
    exponente = (pp - 1) * (qq - 1) // gcd(pp - 1, qq - 1)   # lcm(p'-1, q'-1)
    return _orden_multiplicativo(2, n, exponente), lam
