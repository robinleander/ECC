import random
import gmpy2
from gmpy2 import mpz

class NotPrime(Exception):
    def __init__(self, msg):
        super().__init__(msg)


class NotOnCurve(Exception):
    def __init__(self, msg):
        super().__init__(msg)


class Singular(Exception):
    def __init__(self, msg):
        super().__init__(msg)


def is_quadratic_residue(n, p):
    """
    Returns true iff $\exists s*s = n \in Z/pZ$ where p is prime
    """
    return gmpy2.powmod(n, (p-1)//2, p) == 1


def isqrt_modp(n, p):
    """
    Compute an integer s such that s*s = n in Z/pZ
    or return None if there is no such integer
    Implements Tonelli-Shanks algorithm

    ```
    assert isqrt_modp(5, 41) == 28
    p = None
    while p is None:
        p = rand_prime(3, 1000000)
    for x in range(p):
        s = isqrt_modp(x, p)
        if s is not None:
            assert gmpy2.powmod(s, 2, p) == x
```
    """
    if not is_quadratic_residue(n, p):
        return None
    Q, S = gmpy2.remove(p-1, 2)
    for z in range(3, p):
        if not is_quadratic_residue(z, p):
            break
    M = S
    c = gmpy2.powmod(z, Q, p)
    t = gmpy2.powmod(n, Q, p)
    R = gmpy2.powmod(n, (Q+1)//2, p)
    for its in range(10000):
        if t == 0:
            return 0
        if t == 1:
            return R
        for i in range(1, M):
            if t == 1:
                break
            t = gmpy2.powmod(t, 2, p)
        b = gmpy2.powmod(c, gmpy2.powmod(2, M-i-1, p), p)
        M = i
        b2 = gmpy2.powmod(b, 2, p)
        c = b2
        t = gmpy2.f_mod(t * b2, p)
        R = gmpy2.f_mod(R * b, p)


def rand_prime(start, stop):
    """
    Generates a random prime $\in {start, ..., stop-1}$
    """
    p = random.randrange(start, stop)
    for i in range(p, stop):
        if gmpy2.is_prime(i):
            return i


class Zp:
    """
    Z_p
    """
    def __init__(self, p=None):
        if p is None:
            p = rand_prime(3, 2**256)
        p = mpz(p)
        if not gmpy2.is_prime(p):
            raise NotPrime('{} is not prime'.format(p))
        self.p = p

    def __call__(self, x):
        return Value(self, x)

    def __repr__(self):
        return '{}({})'.format(type(self).__name__, self.p)

    def __eq__(self, b):
        return type(b) == type(self)  and self.p == b.p

    def random_element(self):
        return self(random.randrange(0, self.p))

class Value:
    def __init__(self, zp, v):
        if type(zp) != Zp:
            zp = Zp(zp)
        assert type(zp) == Zp
        self.zp = zp
        if type(v) == Value:
            v = v.v
        self.v = gmpy2.f_mod(v, zp.p)

    def __add__(self, b):
        if type(b) != Value:
            b = Value(self.zp, b)
        assert self.zp == b.zp
        return Value(self.zp, self.v + b.v)

    def __sub__(self, b):
        if type(b) != Value:
            b = Value(self.zp, b)
        assert self.zp == b.zp
        return Value(self.zp, self.v - b.v)

    def __neg__(self):
        return Value(self.zp, -self.v)

    def __mul__(self, b):
        if type(b) != Value:
            b = Value(self.zp, b)
        assert self.zp == b.zp
        return Value(self.zp, self.v * b.v)

    def __floordiv__(self, b):
        if type(b) != Value:
            b = Value(self.zp, b)
        assert self.zp == b.zp
        return Value(self.zp, self * b.inverse())

    def __eq__(self, b):
        if type(b) != Value:
            b = Value(self.zp, b)
        return type(b) == type(self) and self.zp == b.zp and self.v == b.v

    def inverse(self):
        return Value(self.zp, gmpy2.invert(self.v, self.zp.p))

    def __repr__(self):
        return '{}({}, {})'.format(type(self).__name__, repr(self.zp), self.v)

    def __str__(self):
        return str(self.v)


class EC:
    """
    y^2 = x^3 + a*x + b
    """
    def __init__(self, zp, a, b):
        if type(zp) != Zp:
            zp = Zp(zp)
        assert type(zp) == Zp
        if type(a) != Value:
            a = Value(zp, a)
        else:
            assert a.zp == zp
        if type(b) != Value:
            b = Value(zp, b)
        else:
            assert b.zp == zp
        if a*a*a*4 + b*b*27 == 0:
            raise Singular('{} is singular'.format(self))
        self.zp = zp
        self.a = a
        self.b = b

    def at(self, x):
        x = Value(self.zp, x)
        return x * x * x + self.a * x + self.b

    def contains(self, x, y):
        y = Value(self.zp, y)
        return y * y == self.at(x)

    def __call__(self, x, y):
        return P(self, x, y)

    def __repr__(self):
        return '{}({}, {}, {})'.format(
                type(self).__name__,
                repr(self.zp),
                repr(self.a),
                repr(self.b)
                )

    def __str__(self):
        return '{}({}, {}, {})'.format(type(self).__name__, self.zp.p, self.a, self.b)

    def __eq__(self, b):
        return type(b) == type(self) and self.zp == b.zp and self.a == b.a and self.b == b.b

    """
    Generate a random point on the curve
    """
    def random_element(self):
        for its in range(10000):
            try:
                x = self.zp.random_element()
                y2 = self.at(x)
                y = isqrt_modp(y2.v, self.zp.p)
                if y is None:
                    continue
                return P(self, x, random.choice([-1, 1]) * y)
            except NotOnCurve as e:
                pass


class Point:
    def __init__(self, E):
        assert type(E) == EC
        self.E = E


class P(Point):
    def __init__(self, E, x, y):
        super().__init__(E)
        self.x = Value(E.zp, x)
        self.y = Value(E.zp, y)
        if not E.contains(x, y):
            raise NotOnCurve('{} is not on the curve'.format(self))

    def __repr__(self):
        return '{}({}, {}, {})'.format(type(self).__name__, repr(self.E), repr(self.x), repr(self.y))
    
    def __str__(self):
        return '{}({}, {})'.format(type(self).__name__, self.x, self.y)

    def __eq__(self, b):
        return type(b) == type(self) and self.E == b.E and self.x == b.x and self.y == b.y

    def __add__(self, b):
        if type(b) == O:
            return b
        if self.x == b.x and self.y == -b.y:
            return O(self.E)
        if self == b:
            s = (self.x*self.x*3 + self.E.a)//(self.y*2)
        else:
            s = (b.y - self.y)//(b.x - self.x)
        x = s*s - self.x - b.x
        y = s*(self.x - x) - self.y
        return P(self.E, x, y)

    def __sub__(self, b):
        return self + (-b)

    def __neg__(self):
        return P(self.E, self.x, -self.y)

    def __mul__(self, b):
        assert type(b) != P and type(b) != O
        rv = O(self.E)
        fm = self
        b = mpz(b)
        for i in range(b.bit_length()):
            if b.bit_test(i):
                rv = rv + fm
            fm += fm
        return rv


class O(Point):
    """
    The neutral element of the EC group
    """
    def __init__(self, E):
        super().__init__(E)

    def __repr__(self):
        return '{}({})'.format(type(self).__name__, self.E)

    def __str__(self):
        return '{}()'.format(type(self).__name__)

    def __eq__(self, b):
        return type(b) == type(self) and self.E == b.E

    def __add__(self, b):
        return b

    def __mul__(self, b):
        return b


brainpoolP256r1 = EC(0xA9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5377, 0x7D5A0975FC2C3057EEF67530417AFFE7FB8055C126DC5C6CE94A4B44F330B5D9, 0x26DC5C6CE94A4B44F330B5D9BBD77CBF958416295CF7E1CE6BCCDC18FF8C07B6)
brainpoolP256r1_g = brainpoolP256r1(0x8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262, 0x547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997)
brainpoolP256r1_n = 0xA9FB57DBA1EEA9BC3E660A909D838D718C397AA3B561A6F7901E0E82974856A7
