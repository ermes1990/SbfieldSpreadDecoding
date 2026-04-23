#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Field F_64 = GF(2^6)
F.<a> = GF(2^6)

# Field F_8
K = GF(2^3, name='b')

# Embedding of F_8 in F_64
# Let’s take an element of order 7 (a multiplicative generator of F_8^*)
b = a^9  # a^(63/7) = a^9 has order 7

F8_elements = [0] + [b^i for i in range(7)]


GREEN = "\033[92m"
RESET = "\033[0m"

print("Elements of F_64 (in green if elements of F_8):")

for x in F:
    if x in F8_elements:
        print(f"{GREEN}{x}{RESET}")
    else:
        print(x)

F.modulus()

V = VectorSpace(GF(2), 6)


gens = [V((F(1)).list()), V((a^3).list()), V((a^4).list()), V((a^5).list())]

calV = V.subspace(gens)

calV_elements = [F(v.list()) for v in calV.list()]


GREEN = "\033[92m"
RESET = "\033[0m"

print("\n\nElements of \calV (if green, they belong to F_8):")
for x in calV_elements:
    if x in F8_elements:
        print(f"{GREEN}{x}{RESET}")
    else:
        print(x)
# x_1calV
x_1calV = calV_elements

m = a^3 + a^2 + a
# x_2calV
x_2calV = [x * m for x in calV_elements]


m = a^3 + a^4
# x_3calV
x_3calV = [x * m for x in calV_elements]

# x_1calV\capx_2calV\capx_3calV
intersection = list(set(x_1calV) & set(x_2calV) & set(x_3calV))

# print x_1calV\capx_2calV\capx_3calV
print("\n \n x_1calV\capx_2calV\capx_3calV:", intersection)

