"""Figure for the note on Banks--Freiberg--Maynard (Proc. LMS 113 (2016); arXiv:1404.5094).

Run:  python make_figures.py   ->  ../figures/fig_bfm_2016_audit.pdf

(a) The leading part f(x) = x - m - 2 w x^2 / B of the normalised main term of S' in Theorem 4.3(ii),
    x = rho*delta*log k, for m = 1, 2, 3: the paper's (w, B) = (m+1, 8m^2+8m) (touches 0 at x = 2m),
    (m+1, 8m^2+8m+1), and the repair (m, 8m^2+1) of the note.  The O(delta) and O((log k)^{-1/2})
    terms are not explicit in the paper and are not plotted.
(b) The bracket of (5.16) in Lemma 5.2 as a function of c, where 2z = c y log_2 y / log_3 y:
    the limit 1 - 4 c e^{-gamma} and the finite-y values along y = exp exp exp t (formula as in
    verify_bfm_2016.py).
"""
import os
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mp.mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "fig_bfm_2016_audit.pdf")

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 8.5,
    "axes.labelsize": 9, "legend.fontsize": 7.6, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 0.7, "pdf.fonttype": 42,
})
# Okabe-Ito
COL = {1: "#0072B2", 2: "#D55E00", 3: "#009E73"}
GREY = "#555555"

fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.3, 2.75), gridspec_kw={"wspace": 0.30})

# ---------------------------------------------------------------- (a)
for m in (1, 2, 3):
    x = np.linspace(2 * m - 1.35, 2 * m + 1.35, 400)
    variants = [((m + 1, 8 * m * m + 8 * m), "-", 1.3),
                ((m + 1, 8 * m * m + 8 * m + 1), "--", 1.0),
                ((m, 8 * m * m + 1), ":", 1.5)]
    for (w, B), ls, lw in variants:
        ax.plot(x, x - m - 2 * w * x ** 2 / B, ls=ls, lw=lw, color=COL[m])
    ax.plot([2 * m], [0], "o", ms=3.2, color=COL[m], zorder=5)
    ax.text(2 * m, -0.47, rf"$m={m}$", ha="center", va="bottom", color=COL[m], fontsize=8)
ax.axhline(0, color=GREY, lw=0.6)
ax.set_xlim(0.2, 7.8)
ax.set_ylim(-0.5, 0.36)
ax.set_xlabel(r"$x=\rho\delta\log k$")
ax.set_ylabel(r"$x-m-2wx^2/B$")
from matplotlib.lines import Line2D
h = [Line2D([], [], color="k", ls="-", lw=1.3), Line2D([], [], color="k", ls="--", lw=1.0),
     Line2D([], [], color="k", ls=":", lw=1.5)]
ax.legend(h, [r"$w=m+1,\ B=8m^2+8m$", r"$w=m+1,\ B=8m^2+8m+1$", r"$w=m,\ B=8m^2+1$"],
          loc="upper right", frameon=False, handlelength=2.2, borderaxespad=0.2)
ax.set_title("(a)", fontsize=9, loc="left")

# ---------------------------------------------------------------- (b)
G = mp.euler
eg = float(mp.exp(-G))


def bracket(t, c):
    L2 = mp.exp(t)
    zy = c * L2 / (2 * t)
    return float(1 - 2 * mp.exp(-G) * zy * mp.log(zy * t) / (L2 / 4))


c = np.linspace(0.02, 1, 200)
for t, col, ls in ((5, "#CC79A7", "--"), (12, "#E69F00", "-.")):
    bx.plot(c, [bracket(t, mp.mpf(ci)) for ci in c], ls=ls, lw=1.0, color=col,
            label=rf"$y=\exp\exp\exp {t}$")
bx.plot(c, 1 - 4 * c * eg, color="#0072B2", lw=1.5, label=r"limit $1-4ce^{-\gamma}$")
bx.axhline(0, color=GREY, lw=0.6)
bx.axhline(1 - eg, color=GREY, lw=0.7, ls=":")
cz = float(mp.exp(G) / 4)
bx.plot([cz], [0], "o", ms=3.5, color="#0072B2")
bx.annotate(r"$c=e^{\gamma}/4$", xy=(cz, 0), xytext=(0.13, -0.38), fontsize=8,
            arrowprops=dict(arrowstyle="-", lw=0.6, color=GREY))
bx.plot([0.25], [1 - eg], "s", ms=3.5, color="k")
bx.text(0.985, 1 - eg + 0.04, r"printed $1-e^{-\gamma}$ (level)", fontsize=8, va="bottom", ha="right")
bx.plot([1], [1 - 4 * eg], "D", ms=3.5, color="#D55E00", clip_on=False, zorder=5)
bx.annotate(r"$\S6$: $c=1$", xy=(1, 1 - 4 * eg), xytext=(0.66, -1.12), fontsize=8,
            arrowprops=dict(arrowstyle="-", lw=0.6, color=GREY))
bx.set_xlim(0, 1)
bx.set_ylim(-1.4, 1.05)
bx.set_xlabel(r"$c=2z\log_3 y\,/\,(y\log_2 y)$")
bx.set_ylabel("bracket in (5.16)")
bx.legend(loc="lower left", frameon=False, borderaxespad=0.2)
bx.set_title("(b)", fontsize=9, loc="left")

fig.subplots_adjust(left=0.095, right=0.985, bottom=0.16, top=0.92)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT)
print("wrote", OUT)



