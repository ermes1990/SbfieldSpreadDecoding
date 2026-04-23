#!/usr/bin/env python
# coding: utf-8

# In[1]:


import random as pyrandom
from sage.all import *

from IPython.display import display, HTML

def col(x, color):
    return f"<span style='color:{color}'>{x}</span>"

print_block = """
    <table style="border-collapse:collapse; font-family:monospace;">
    <tr>
    <th style="color:cyan">n</th>
    <th style="color:orange">k</th>
    <th style="color:green">i</th>
    <th style="color:magenta">t</th>

    <th style="color:blue">zr≤LB_GR</th>
    <th style="color:red">zr_SLT_OurB</th>
    <th style="color:cyan">LB_GR=OurB</th>
    <th style="color:orange">OurB_SLT_UP_GR</th>
    <th style="color:green">OurB_SLT_tot</th>
    <th style="color:magenta">UP_GR_SLT_tot</th>
    <th style="color:cyan">(i-t)≤(n-k)</th>
    </tr>
"""

print_block1 = """
    <table style="border-collapse:collapse; font-family:monospace;">
    <tr>
    <th style="color:cyan">n</th>
    <th style="color:orange">k</th>
    <th style="color:green">i</th>
    <th style="color:magenta">t</th>

    <th style="color:blue">zr_SLT_LB_GR==OurB==UP_GR_SLT_tot</th>
    <th style="color:blue">tot-bounds</th>
    <th style="color:blue">(tot-bounds)/tot</th>
    </tr>
"""

print_block2 = """
    <table style="border-collapse:collapse; font-family:monospace;">
    <tr>
    <th style="color:cyan">n</th>
    <th style="color:orange">k</th>
    <th style="color:green">i</th>
    <th style="color:magenta">t</th>

    <th style="color:blue">denominator</th>
    <th style="color:red">delta</th>
    <th style="color:cyan">omega</th>
    <th style="color:orange">zr≤OurB</th>
    <th style="color:green">OurB==tot</th>
    <th style="color:cyan">(i-t)≤(n-k)</th>
    </tr>
"""

ct = 0
rg = 100

for _ in range(rg):
    # parameters q=p^e
    p = pyrandom.choice(list(primes(2,10)))
    e = pyrandom.randint(10,20)
    q = p^e

    # parameters n,k,i,t
    n = pyrandom.randint(20,100)
    while is_odd(n):
        n = pyrandom.randint(20,100)

    divisors = [d for d in range(2, n) if n % d == 0]
    k = pyrandom.choice(divisors)

    t = pyrandom.randint(1, k-1)
    b = min (floor((t*n-k*t+t^2+t)/(t+1)), floor((t*n)/(t+1)))
    c = floor(n/4)
    i = pyrandom.randint(1, b)


    sum1 = sum(
        q^((k-h)*(i-h)) *
        q_binomial(k, h, q) *
        q_binomial(n-k, i-h, q)
        for h in range(t+1)
    )

    our_bound = q_binomial(n, i, q) - (((q^n - 1)/(q^k - 1)) * (q_binomial(n, i, q) - sum1))

    total_subspaces_fixed_dimension_i = q_binomial(n, i, q)

    # ----------------------------
    # delta_q
    # ----------------------------
    delta = sum(
        sum(
            q_binomial(k, l, q)
            * q_binomial(k-l, b2-l, q)
            * q_binomial(n-b2, i-b2, q)
            * (-1)^(b2-l)
            * q^(binomial(b2-l, 2))
            for b2 in range(l, k+1)
        )
        for l in range(t+1, k+1)
    )

    # ----------------------------
    # omega_q
    # ----------------------------
    omega = sum(
        q_binomial(k, l, q)
        * q_binomial(k, lp, q)
        * sum(
            sum(
                q_binomial(k-l, r-l, q)
                * q_binomial(k-lp, s-lp, q)
                * q_binomial(n-r-s, i-r-s, q)
                * (-1)^(r+s-l-lp)
                * q^(binomial(r-l, 2)+binomial(s-lp, 2))
                for s in range(lp, k+1)
            )
            for r in range(l, k+1)
        )
        for lp in range(t+1, k+1)
        for l in range(t+1, k+1)
    )

    # ----------------------------
    # bound at least & at most
    # ----------------------------
    lower_bound_GR = q_binomial(n, i, q) - ((q^n - 1)/(q^k - 1)) * delta

    denominator = (delta + ((((q^n - 1)/(q^k - 1)) - 1) * omega))
    if denominator > 0:
        bound_at_most_int = (
            (((q^n - 1)/(q^k - 1)) * delta^2)
            / denominator
        )

        upper_bound_GR = q_binomial(n, i, q) - bound_at_most_int

        cond0 = 0 <= lower_bound_GR
        ourboundSTRICTLYGRzero = 0 < our_bound
        cond1 = lower_bound_GR == our_bound
        cond2 = our_bound < upper_bound_GR
        cond3 = upper_bound_GR < total_subspaces_fixed_dimension_i
        cond4 = our_bound < total_subspaces_fixed_dimension_i
        cond5 = (i-t) <= (n-k)
        cond6 = (0 < our_bound) and ( our_bound == lower_bound_GR) and (our_bound == upper_bound_GR) and (our_bound < total_subspaces_fixed_dimension_i)
        diff = total_subspaces_fixed_dimension_i - upper_bound_GR
        diffOVERtot = diff /total_subspaces_fixed_dimension_i

        if cond2: 
            row = f"""
                <tr>
                <td align="right">{col(n,'cyan')}</td>
                <td align="right">{col(k,'orange')}</td>
                <td align="right">{col(i,'green')}</td>
                <td align="right">{col(t,'magenta')}</td>

                <td align="center">{col('✓' if cond0 else '!', 'blue')}</td>
                <td align="center">{col('✓' if ourboundSTRICTLYGRzero else '!', 'red')}</td>
                <td align="center">{col('✓' if cond1 else '!', 'cyan')}</td>
                <td align="center">{col('✓' if cond2 else '!', 'orange')}</td>
                <td align="center">{col('✓' if cond4 else '!', 'green')}</td>
                <td align="center">{col('✓' if cond3 else '!', 'magenta')}</td>
                <td align="center">{col('✓' if cond5 else '!', 'cyan')}</td>
                </tr>
            """
            print_block += row
        else:
            row1 = f"""
                <tr>
                <td align="right">{col(n,'cyan')}</td>
                <td align="right">{col(k,'orange')}</td>
                <td align="right">{col(i,'green')}</td>
                <td align="right">{col(t,'magenta')}</td>

                <td align="center">{col('✓' if cond6 else '!', 'blue')}</td>
                <td align="center">{col(diff.n(digits = 5), 'blue')}</td>
                <td align="center">{col(diffOVERtot.n(digits = 5), 'blue')}</td>
                </tr>
            """
            print_block1 += row1


    else:
        ourboundSTRICTLYGRzero = 0 < our_bound
        cond4 = our_bound == total_subspaces_fixed_dimension_i
        cond5 = (i-t) <= (n-k)
        row2 = f"""
            <tr>
            <td align="right">{col(n,'cyan')}</td>
            <td align="right">{col(k,'orange')}</td>
            <td align="right">{col(i,'green')}</td>
            <td align="right">{col(t,'magenta')}</td>

            <td align="center">{col(denominator, 'blue')}</td>
            <td align="center">{col(delta, 'red')}</td>
            <td align="center">{col(omega, 'cyan')}</td>
            <td align="center">{col('✓' if ourboundSTRICTLYGRzero else '!', 'orange')}</td>
            <td align="center">{col('✓' if cond4 else '!', 'green')}</td>
            <td align="center">{col('✓' if cond5 else '!', 'cyan')}</td>
            </tr>
        """
        print_block2 += row2

display(HTML(print_block + "</table>"))
print("\n Exact estimate")
display(HTML(print_block1 + "</table>"))
print("\n DENOMINATOR = 0")
display(HTML(print_block2 + "</table>"))

