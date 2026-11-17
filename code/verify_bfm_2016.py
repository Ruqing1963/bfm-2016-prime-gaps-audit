"""Verification for the note on Banks--Freiberg--Maynard, "On limit points of the sequence of
normalized prime gaps" (Proc. LMS 113 (2016); arXiv:1404.5094).  Line numbers refer to the arXiv source.

Run:  python verify_bfm_2016.py   ->  ../data/verify_bfm_2016.txt

Sections
  1. Theorem 4.3(i): normalised main term with unordered / ordered pairs (exact rationals).
  2. Theorem 4.3(ii): main term x - m - 2 w x^2 / B, x = rho*delta*log k; maxima and the tie at
     (w, B) = (m+1, 8m^2+8m); the values for B = 8m^2+8m+1 and for the repair (w, B) = (m, 8m^2+1).
  3. Exhaustive check of the combinatorial step (l.1931-1942): for every sequence of block counts
     c_1..c_t in {1,2,3} (zero blocks play no role), does a positive weight
        sum c_i - m - w * sum C(c_i, 2) > 0
     force (P) m+1 consecutive non-empty blocks that are all singletons (the paper's (4.21)), or
     (Q) m+1 consecutive non-empty blocks whose m-1 inner blocks are singletons (what Section 6 needs)?
  4. Lemma 5.2: the bracket 1 - 2 e^{-gamma} z log(z/y_2) / (y log y_1), y_1 = (log y)^{1/4},
     y_2 = y / log_3 y, 2z = c y log_2 y / log_3 y, along y = exp exp exp t, and its limit 1 - 4 c e^{-gamma}.
  5. Corollary 1.2 constants and the k-multiple at l.2688.
"""
import itertools
import os
from fractions import Fraction as Fr
from math import comb

import mpmath as mp

mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "verify_bfm_2016.txt")
lines = []


def out(s=""):
    print(s)
    lines.append(s)


# ------------------------------------------------------------------ 1. part (i)
out("1. Theorem 4.3(i): x = 2m, 8m+1 blocks of size k/(8m+1), pair constant 4, weight 1")
out("   main term = x - m - 4 (x/k)^2 * #pairs ;  unordered: #pairs = k(k/(8m+1)-1)/2")
for m in range(1, 6):
    for k in (Fr((8 * m + 1) * 100), Fr((8 * m + 1) * 10 ** 6)):
        x = Fr(2 * m)
        nb = 8 * m + 1
        unord = k * (k / nb - 1) / 2
        mt_u = x - m - 4 * (x / k) ** 2 * unord
        mt_o = x - m - 4 * (x / k) ** 2 * (2 * unord)
        claim = Fr(m, 8 * m + 1) + Fr(8 * m * m) / k
        out(f"   m={m} k={int(k):>9}: unordered {float(mt_u):+.10f} (paper m/(8m+1)+8m^2/k = {float(claim):+.10f},"
            f" equal: {mt_u == claim});  ordered {float(mt_o):+.6f} = (m-8m^2)/(8m+1)+16m^2/k: "
            f"{mt_o == Fr(m - 8 * m * m, 8 * m + 1) + Fr(16 * m * m) / k}")
out()

# ------------------------------------------------------------------ 2. part (ii)
out("2. Theorem 4.3(ii): f(x) = x - m - 2 w x^2 / B  (leading part, k -> infinity)")
out("   max_x f = B/(8w) - m at x = B/(4w);  positive iff B > 8 m w")
for m in range(1, 7):
    rows = []
    for (w, B, tag) in ((m + 1, 8 * m * m + 8 * m, "paper  w=m+1,B=8m^2+8m  "),
                        (m + 1, 8 * m * m + 8 * m + 1, "        w=m+1,B=8m^2+8m+1"),
                        (m + 2, 8 * m * m + 16 * m + 1, "        w=m+2,B=8m^2+16m+1"),
                        (m, 8 * m * m + 1, "repair w=m,  B=8m^2+1   ")):
        w, B = Fr(w), Fr(B)
        f2m = 2 * m - m - 2 * w * (2 * m) ** 2 / B
        fmax = B / (8 * w) - m
        out(f"   m={m} {tag}: f(2m) = {str(f2m):>10} ;  max = {str(fmax):>10} at x = {str(B / (4 * w))}")
    out(f"        leftover 2(m+1)(2m)^2/k at the tie = {8 * m * m * (m + 1)}/k")
out()

# ------------------------------------------------------------------ 3. combinatorics


def has_P(c, m):
    """m+1 consecutive non-empty blocks all equal to 1."""
    run = 0
    for ci in c:
        run = run + 1 if ci == 1 else 0
        if run >= m + 1:
            return True
    return False


def has_Q(c, m):
    """window of m+1 consecutive non-empty blocks with inner m-1 blocks equal to 1."""
    t = len(c)
    for s in range(0, t - m):
        if all(ci == 1 for ci in c[s + 1:s + m]):
            return True
    return False


out("3. Combinatorial step: sequences of non-empty block counts c_i in {1,2,3}, length t <= TMAX")
TMAX = 13
for m in (2, 3, 4):
    for w in (Fr(m) - Fr(1, 2), Fr(m), Fr(m + 1), Fr(m + 2)):
        nP = nQ = npos = 0
        exP = exQ = None
        for t in range(1, TMAX + 1):
            for c in itertools.product((1, 2, 3), repeat=t):
                val = sum(c) - m - w * sum(comb(ci, 2) for ci in c)
                if val <= 0:
                    continue
                npos += 1
                if not has_P(c, m):
                    nP += 1
                    if exP is None:
                        exP = c
                if not has_Q(c, m):
                    nQ += 1
                    if exQ is None:
                        exQ = c
        out(f"   m={m} w={str(w):>4}: positive sequences {npos:>8}; (4.21) fails {nP:>8}"
            f" (first {exP}); window fails {nQ:>6} (first {exQ})")
m = 2
c = (1, 1, 2, 1, 1)
out(f"   example m=2, w=m+1=3, blocks S S D S S = {c}: weight = {sum(c)} - {m} - 3*1 = "
    f"{sum(c) - m - 3 * sum(comb(ci, 2) for ci in c)} > 0; (4.21): {has_P(c, m)}; window: {has_Q(c, m)}")
out("   paper's count (l.1933): #singletons >= m+1+m r ; correct for w=m+1: >= m+1+(m-1) r")
out()

# ------------------------------------------------------------------ 4. bracket
G = mp.euler
eg = mp.exp(-G)


def bracket(t, c):
    L2 = mp.exp(t)            # log_2 y
    L3 = mp.mpf(t)            # log_3 y
    zy = c * L2 / (2 * L3)    # z / y
    return 1 - 2 * eg * zy * mp.log(zy * L3) / (L2 / 4)


out("4. Lemma 5.2 bracket, 2z = c y log_2 y / log_3 y, y = exp exp exp t")
out(f"   e^-gamma = {mp.nstr(eg, 12)};  1 - 4e^-gamma = {mp.nstr(1 - 4 * eg, 12)};"
    f"  zero c* = e^gamma/4 = {mp.nstr(mp.exp(G) / 4, 12)};  printed 1 - e^-gamma = {mp.nstr(1 - eg, 12)}")
for c in (mp.mpf(1) / 4, mp.exp(G) / 4, mp.mpf(1) / 2, mp.mpf(1)):
    vals = ", ".join(f"t={t}: {mp.nstr(bracket(t, c), 6)}" for t in (4, 8, 20, 80, 400))
    out(f"   c={mp.nstr(c, 6):>8}: limit {mp.nstr(1 - 4 * c * eg, 8):>11};  {vals}")
out("   (z/y_2 = c log_2 y / 2, so log(z/y_2) = log_3 y + log(c/2); log y_1 = log_2 y / 4)")
# conditions used at l.2338 and (5.10) at the top of the range
t = 8
L2 = mp.exp(t)
out(f"   side conditions at t={t}, c=1: log(z/y_2) = {mp.nstr(mp.log(L2 / 2), 6)} < log y_1 = {mp.nstr(L2 / 4, 6)}")
out()

# ------------------------------------------------------------------ 5. misc
H8 = sum(Fr(1, j) for j in range(1, 9))
out("5. Corollary 1.2: kappa - 1 <= 8 -> T/8 = 12.5%;  effective constant 1/(8 H_8) = "
    f"{float(1 / (8 * H8)):.6f} > 1/22 = {1 / 22:.6f}")
for m in (1,):
    req = (8 * m + 1) * (8 * m * m + 8 * m)
    out(f"   l.2688: 9*17 = 153; (8m+1)(8m^2+8m) = {req} for m=1; 153 % {req} = {153 % req};"
        f" 153 % 9 = {153 % 9}")
out("   repaired Theorem 4.3: k a multiple of (8m+1)(8m^2+1): " +
    ", ".join(f"m={m}: {(8 * m + 1) * (8 * m * m + 1)}" for m in range(1, 5)))

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
