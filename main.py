#!/usr/bin/env python3

from ecc.ec import *
import random

E = EC(11, 3, 7)

g = brainpoolP256r1_g
a = random.randrange(2, 2**256)
b = random.randrange(2, 2**256)

print("Alice chooses secret a = {}".format(a))
print("  Bob chooses secret b = {}".format(b))
print()

A = brainpoolP256r1_g * a
B = brainpoolP256r1_g * b
print("Alice sends Bob A = a*G = {}".format(A))
print("Bob sends Alice B = b*G = {}".format(B))
print()

Sa = B * a
Sb = A * b
print("Both parties compute the shared secret")
print("Alice computes S = a*B = a*b*G = {}".format(Sa))
print("  Bob computes S = b*A = b*a*G = {}".format(Sb))
print()

if Sa == Sb:
    print("The secret is shared.")
else:
    raise Exception("The secret is different!!")

