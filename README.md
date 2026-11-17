# A note on limit points of normalized prime gaps

Ruqing Chen (GUT Geoservice Inc., Montreal). A 4-page note on W. D. Banks, T. Freiberg, J. Maynard,
*On limit points of the sequence of normalized prime gaps*, Proc. London Math. Soc. (3) 113 (2016), no. 4,
515–539, DOI 10.1112/plms/pdw036, arXiv:1404.5094. Line numbers refer to the
arXiv source of that paper (https://arxiv.org/abs/1404.5094, "Download source"), which is not redistributed here.

Archived at Zenodo: DOI [10.5281/zenodo.24315678](https://doi.org/10.5281/zenodo.24315678) (reserved).

## Files
- `bfm_prime_gaps_audit.tex`, `bfm_prime_gaps_audit.pdf`: the note. pdflatex compiles it with no warnings.
- `code/verify_bfm_2016.py`: the verification script (fractions + mpmath, a few minutes because of the exhaustive check).
- `data/verify_bfm_2016.txt`: the script's output.
- `code/make_figures.py` -> `figures/fig_bfm_2016_audit.pdf`.
## Findings
1. **Theorem 4.3(ii) has a tie.** With weight m+1, B = 8m²+8m blocks and x = ρδ log k = 2m, the main term is
   x − m − x²/(4m) = −(x−2m)²/(4m), so its maximum is exactly 0. The leftover 8m²(m+1)/k is smaller than the O(δ) error.
   Part (i) is correct: its main term m/(8m+1) + 8m²/k uses unordered pairs.
2. **Its combinatorial step is wrong for weight m+1.** At l.1933 the paper claims at least m+1+mr singleton blocks, but
   with weight m+1 the bound is only m+1+(m−1)r. Counterexample (m=2): blocks (1,1,2,1,1). So changing B to
   8m²+8m+1 alone does not fix the proof. Weight m+2 does give (4.21), but needs B ≥ 8m²+16m+1. This is Lau's
   threshold at r = 2.
3. **Repair (Lemma 1, Proposition 2 of the note).** Section 6 needs only m+1 *consecutive* primes in distinct blocks.
   That follows from weight m, and with weight m, B = 8m²+1 suffices. So Theorem 1.3 holds as stated, and in fact for
   any 8m²+1 reals. This needs the h_i to be chosen increasing in i, so that the blocks are ordered as integers.
4. **l.2688:** 9×17 is not a multiple of 9×16. This is harmless, because only part (i) is used (9 | k). Also, (4.20)
   should say ≥ 1, not = 1.
5. **Lemma 5.2 bracket.** Write 2z = c·y log₂y/log₃y. The bracket then tends to 1 − 4c e^{−γ}. At c = 1, the value used
   in §6, the limit is −1.2458. It is negative for c > e^γ/4 = 0.4453. The printed 1 − e^{−γ} is the value at c = 1/4.
   Repair: require 2z ≤ c₀ y log₂y/log₃y with c₀ < e^γ/4, and in §6 take z ≍ y. The remark at l.2744 survives up to a
   constant factor. Merikoski (JLMS 2020) quotes the lemma verbatim with c = 1, and the same repair applies to it.
6. Misprints: l.1813, l.2389, l.2645, l.2697. The flagged l.287 (j ≤ κ) is not an error; only the hat notation is off.
