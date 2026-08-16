"""Generate THEORY.pdf: the physics, the error budget, and what each setup can do.

    python make_theory_pdf.py

Development script; the only things in this repository needing matplotlib are this and
make_field_card.py. The field package itself stays pure standard library.

Five pages:
  1  the measurement problem and the governing equations
  2  WHAT YOU HAVE decides WHAT YOU GET  (the capability chart)
  3  the error budget, quantified
  4  two instruments, and how to cross-check them
  5  derived products the absolute scale supports
"""

import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Arc, Circle, FancyArrowPatch, FancyBboxPatch, Polygon

OUT = os.path.dirname(os.path.abspath(__file__))

SUN, WATER, SKY = "#f5a623", "#1f7a99", "#7fb3d5"
PANEL, GOOD, BAD, INK = "#d9534f", "#2e7d32", "#c0392b", "#1a1a1a"
ED = "#6a3d9a"
PAGE = (11.7, 16.5)          # A3 portrait


def overlay(fig):
    """A full-page axes in 0-1 coordinates.

    set_xlim/set_ylim are ESSENTIAL: any ax.plot on this axes would otherwise
    autoscale the data limits and silently move every text placement on the page.
    """
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return ax


#: one 9 pt line at linespacing 1.42, as a fraction of the 16.5 in page height
LINE = 12.8 / (16.5 * 72.0)


class Flow:
    """Lay out text blocks top-down inside a box.

    Placing an equation at a hand-chosen y and the prose at another is how equations
    end up printed THROUGH body text; this advances one cursor instead.
    """

    def __init__(self, ax, x, top):
        self.ax, self.x, self.y = ax, x, top

    def text(self, body, size=9.0, color=INK, dy=None):
        self.ax.text(self.x, self.y, body, fontsize=size, color=color, va="top",
                     ha="left", linespacing=1.42, zorder=3)
        self.y -= (dy if dy is not None else body.count("\n") + 1) * LINE * size / 9.0
        return self

    def eq(self, tex, size=14.0, xc=0.5, pad=0.6):
        self.y -= pad * LINE
        self.ax.text(xc, self.y, tex, fontsize=size, color=INK, va="top", ha="center",
                     zorder=3)
        self.y -= (size / 9.0) * LINE * 1.9 + pad * LINE
        return self

    def gap(self, n=1.0):
        self.y -= n * LINE
        return self


def header(fig, n, title, subtitle):
    fig.text(0.5, 0.965, title, fontsize=19, weight="bold", ha="center", color=INK)
    fig.text(0.5, 0.945, subtitle, fontsize=10.5, ha="center", color="#555")
    fig.text(0.955, 0.012, "page %d of 5" % n, fontsize=8.5, ha="right", color="#888")
    fig.text(0.045, 0.012, "fieldrrs — above-water R$_{rs}$ theory", fontsize=8.5,
             color="#888")


def box(ax, x, y, w, h, title, col, fc, body, tsize=11, bsize=8.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008",
                                fc=fc, ec=col, lw=2.0, transform=ax.transAxes,
                                zorder=2))
    ax.text(x + w / 2, y + h - 0.030, title, fontsize=tsize, weight="bold", color=col,
            ha="center", transform=ax.transAxes, zorder=3)
    ax.text(x + 0.016, y + h - 0.062, body, fontsize=bsize, color=INK, va="top",
            ha="left", transform=ax.transAxes, linespacing=1.42, zorder=3)


# ---------------------------------------------------------------- page 1
def page1():
    fig = plt.figure(figsize=PAGE)
    header(fig, 1, "1 · THE MEASUREMENT PROBLEM",
           "Why three scans, and what each equation is for")
    ax = overlay(fig)

    ax.text(0.055, 0.912,
            "A radiometer pointed at water does NOT measure water-leaving light. It "
            "measures everything\narriving along that line of sight, and the water is "
            "only part of it:",
            fontsize=11, color=INK, va="top")

    ax.text(0.5, 0.866,
            r"$L_t(\lambda)\;=\;L_w(\lambda)\;+\;\rho\,L_{sky}(\lambda)"
            r"\;+\;L_{glint}(\lambda)$",
            fontsize=16, ha="center", color=INK)
    for xf, lab, col in ((0.393, "what you want", GOOD),
                         (0.520, "reflected sky", SKY),
                         (0.638, "specular sun", BAD)):
        ax.annotate("", xy=(xf, 0.860), xytext=(xf, 0.843),
                    arrowprops=dict(arrowstyle="-", lw=1.1, color=col))
        ax.text(xf, 0.839, lab, fontsize=8.4, color=col, ha="center", va="top",
                weight="bold")

    ax.text(0.055, 0.812,
            "The geometry (40° from nadir, 135° in azimuth from the sun) is chosen to "
            "make the third term\nnegligible. The second must be measured and "
            "subtracted, which is why there is a sky scan.",
            fontsize=10, color="#444", va="top")

    # --- the three equations
    ax.add_patch(FancyBboxPatch((0.055, 0.645), 0.89, 0.135,
                                boxstyle="round,pad=0.01", fc="#eef5ff", ec="#8aa8c8",
                                lw=1.6, transform=ax.transAxes))
    ax.text(0.5, 0.760, "THE THREE GOVERNING EQUATIONS", fontsize=11, weight="bold",
            ha="center", color="#33506e")
    ax.text(0.5, 0.729, r"$E_d(\lambda)\;=\;\dfrac{\pi\,L_p(\lambda)}{R_p(\lambda)}$"
                        r"$\qquad\qquad$ (F1)  panel route to irradiance",
            fontsize=13, ha="center", color=INK)
    ax.text(0.5, 0.694, r"$L_w(\lambda)\;=\;L_t(\lambda)\;-\;\rho\,L_{sky}(\lambda)$"
                        r"$\qquad\qquad$ (F2)  remove the reflected sky",
            fontsize=13, ha="center", color=INK)
    ax.text(0.5, 0.660, r"$R_{rs}(\lambda)\;=\;\dfrac{L_w(\lambda)}{E_d(\lambda)}"
                        r"\;-\;\Delta$"
                        r"$\qquad\qquad$ (F3)  normalise, remove residual offset",
            fontsize=13, ha="center", color=INK)

    # --- geometry sketch
    gx, gy, s = 0.30, 0.40, 0.185
    ax.plot([gx - 1.30 * s, gx + 1.15 * s], [gy, gy], color=WATER, lw=2.5)
    ax.fill_between([gx - 1.30 * s, gx + 1.15 * s], gy - 0.09, gy, color=WATER,
                    alpha=0.22)
    ax.plot([gx, gx], [gy, gy + 1.05 * s], color="#bbb", ls=":", lw=1)
    ax.plot([gx - 1.1 * s, gx + 1.1 * s], [gy + 0.5 * s, gy + 0.5 * s], color="#bbb",
            ls="--", lw=1)
    a = math.radians(50.0)
    # length chosen so the WATER ray terminates ON the surface rather than below it
    ray = (0.5 * s) / math.sin(a)
    for sgn, col, lab in ((+1, SKY, "sky  $L_{sky}$"), (-1, WATER, "water  $L_t$")):
        dx, dy = ray * math.cos(a), sgn * ray * math.sin(a)
        ax.annotate("", xy=(gx + dx, gy + 0.5 * s + dy), xytext=(gx, gy + 0.5 * s),
                    arrowprops=dict(arrowstyle="-|>", lw=2.4, color=col))
        ax.text(gx + dx + 0.012, gy + 0.5 * s + dy + sgn * 0.012, lab, fontsize=9.5,
                color=col, weight="bold", va="center")
    ax.text(gx + 0.028, gy + 0.5 * s + 0.052, "40°", fontsize=9, color="#2c6f9b",
            weight="bold")
    ax.text(gx + 0.028, gy + 0.5 * s - 0.062, "40°", fontsize=9, color=WATER,
            weight="bold")
    ax.add_patch(FancyBboxPatch((gx - 0.30 * s, gy + 0.05), 0.16 * s * 2.4, 0.012,
                                boxstyle="round,pad=0.002", fc="white", ec=PANEL, lw=2))
    ax.text(gx - 0.06, gy + 0.075, "panel  $L_p$", fontsize=9.5, color=PANEL,
            weight="bold", ha="center")
    ax.text(gx, gy - 0.135, "SKY and WATER are mirrors through the horizontal:\n"
                            "SAME compass bearing, opposite tilt.",
            fontsize=9, ha="center", color=INK)

    # --- what rho is
    box(ax, 0.545, 0.245, 0.400, 0.300,
        r"$\rho$ — the term everything hinges on", "#8a6000", "#fff6d5",
        "It is NOT the Fresnel coefficient. A flat surface would\n"
        "reflect ~0.021 at 40°. The sea is not flat: it is a\n"
        "distribution of wave facets, each seeing a different\n"
        "part of the sky, so ρ is an EFFECTIVE factor over that\n"
        "distribution.\n\n"
        "That is why it depends on wind — wind sets the facet\n"
        "distribution — and why it depends on the sky, because\n"
        "the facets are sampling the sky.\n\n"
        "ρ = 0.028 (Mobley 1999): 40°/135°, wind < ~5 m/s,\n"
        "clear sky.\n\n"
        "KEY CONSEQUENCE: under a UNIFORM sky there is little\n"
        "to sample between, so facet orientation stops\n"
        "mattering and ρ becomes weakly wind-dependent. This\n"
        "is why overcast is workable.", bsize=8.3)

    ax.text(0.055, 0.205,
            "The residual offset Δ (F3)", fontsize=11.5, weight="bold", color=INK)
    ax.text(0.055, 0.185,
            "Even at the best geometry, (F2) leaves a residual from facets that caught "
            "the sun, from whitecaps,\nand from the sky scan not perfectly matching what "
            "the surface actually reflected. Three treatments:",
            fontsize=9.5, color="#444", va="top")
    ax.text(0.075, 0.140,
            "none                 Δ = 0.  The default, deliberately.\n"
            "nir_zero          assume $R_{rs}$ = 0 over 750–800 nm and subtract the mean.\n"
            "                         VALID in clear oceanic water. In turbid water it "
            "DELETES REAL SIGNAL.\n"
            "nir_similarity   Ruddick et al. (2006): fixed $R_{rs}$(780)/$R_{rs}$(870) = "
            "1.912 solves for Δ.\n"
            "                         Usable in turbid water.",
            fontsize=9.3, color=INK, va="top", linespacing=1.55)
    ax.text(0.055, 0.055,
            "Which one is legal is a property of YOUR WATER, not a software preference. "
            "That is why the default\nis 'none': a correction applied in the field to a "
            "CSV you keep cannot be undone. Keep the raw .sed.",
            fontsize=9.5, color=BAD, va="top")
    return fig


# ---------------------------------------------------------------- page 2
def page2():
    fig = plt.figure(figsize=PAGE)
    header(fig, 2, "2 · WHAT YOU HAVE DECIDES WHAT YOU GET",
           "The exe supports both setups. This is which is which.")
    ax = overlay(fig)

    # ---- SETUP A
    box(ax, 0.045, 0.660, 0.44, 0.250,
        "A ·  NaturaSpec Plus ALONE", WATER, "#eaf4f8",
        "One instrument, narrow field of view (FLENS8 = 8°).\n"
        "THREE radiance scans, one at a time:\n\n"
        "    1. panel      $L_p$        (the REFERENCE scan)\n"
        "    2. sky         $L_{sky}$\n"
        "    3. water     $L_t$\n\n"
        "$E_d$ is INFERRED from the panel by (F1). The panel is\n"
        "a proxy: you never measure $E_d$, you measure $L_p$ and\n"
        "assume the panel converts it faithfully.\n\n"
        "CODE:\n"
        "    rrs_from_sed(water, sky, source='radiance')\n\n"
        "YOU GET:   $R_{rs}$\n"
        "YOU CARRY: panel reflectance, panel levelness,\n"
        "                    panel-to-target TIME LAG, ρ",
        bsize=8.4)

    # ---- SETUP B
    box(ax, 0.515, 0.660, 0.44, 0.250,
        "B ·  + SEPARATE IRRADIANCE SENSOR", ED, "#f2edf9",
        "TWO instruments. The second has a WIDE (hemispherical)\n"
        "field of view and measures $E_d$ directly, logging at the\n"
        "SAME MOMENT as the radiance scans.\n\n"
        "    1. sky         $L_{sky}$      } narrow FOV\n"
        "    2. water     $L_t$          }\n"
        "    3. $E_d$          measured, simultaneous, wide FOV\n\n"
        "No panel in the calculation at all.\n\n"
        "CODE:\n"
        "    rrs_from_separate_ed(wl, Lt, Lsky, Ed, ed_wl)\n\n"
        "YOU GET:   $R_{rs}$, PAR, atmospheric transmittance,\n"
        "                    an $E_d$-stability check\n"
        "YOU CARRY: INTER-INSTRUMENT CALIBRATION, ρ",
        bsize=8.4)

    # arrow between
    ax.annotate("", xy=(0.512, 0.785), xytext=(0.488, 0.785),
                arrowprops=dict(arrowstyle="-|>", lw=2.5, color="#666"))

    # ---- the trade, spelled out
    box(ax, 0.045, 0.470, 0.91, 0.165,
        "THE TRADE — you are swapping one error for another, not removing error",
        INK, "#f7f7f7",
        "A → B  REMOVES three errors:                                    "
        "A → B  ADDS one error:\n"
        "   • the panel-to-target TIME LAG, because $E_d$ is now      "
        "   • $R_{rs}$ divides a radiance from instrument A by an\n"
        "     measured at the same moment as the target.               "
        "     irradiance from instrument B. ANY offset between\n"
        "     Broken cloud stops being disqualifying.                        "
        "     their absolute calibrations is a DIRECT\n"
        "   • the PANEL REFLECTANCE (your 0.99 assumption),         "
        "     MULTIPLICATIVE BIAS on every $R_{rs}$.\n"
        "     a straight multiplicative bias on every number.            "
        "     Measured: a 6 % gain offset → −5.7 % in $R_{rs}$(443),\n"
        "   • PANEL LEVELNESS — the collector defines the              "
        "     uniformly and silently. It does not average out.\n"
        "     horizontal plane itself.",
        tsize=11, bsize=8.3)

    # ---- the cross-check
    box(ax, 0.045, 0.190, 0.91, 0.215,
        "THE CROSS-CHECK — keep the panel, but demote it to a TRANSFER STANDARD",
        GOOD, "#eaf6ea",
        "Point the radiance sensor at the panel while the irradiance sensor sees the "
        "same sky. You now have\nTWO INDEPENDENT ESTIMATES of the same $E_d$, and "
        "their ratio is the inter-instrument factor:\n\n"
        "                    C(λ)  =  [ π · $L_p$(λ) / $R_p$(λ) ]  /  "
        "$E_d^{measured}$(λ)                    (F4)\n\n"
        "    C ≈ 1 and FLAT            the two instruments agree. Use the measured $E_d$ "
        "directly.\n"
        "    C ≠ 1 but FLAT            a gain offset. Decide which absolute scale you "
        "trust before correcting.\n"
        "    C STRUCTURED in λ    NOT a gain. Suspect a stale calibration, a tilted "
        "collector, a shaded panel.\n\n"
        "Run it ONCE PER DEPLOYMENT, and again after anything gets knocked. Reported by\n"
        "cross_calibration_factor(). It is a DIAGNOSTIC first: multiplying C back in "
        "re-introduces the\npanel reflectance that setup B existed to remove.",
        tsize=11, bsize=8.5)
    return fig


# ---------------------------------------------------------------- page 3
def page3():
    fig = plt.figure(figsize=PAGE)
    header(fig, 3, "3 · THE ERROR BUDGET",
           "Which terms actually dominate, with numbers")
    ax = fig.add_axes([0.165, 0.565, 0.79, 0.325])

    terms = ["ρ at 10 m/s\n(clear sky)", "inter-instrument\ncalibration (setup B)",
             "panel reflectance\n(setup A)", "BRDF, uncorrected\n(vs satellite)",
             "panel tilt 5°\n(low sun)", "replicate scatter\n(calm water)"]
    lo = [20, 3, 1, 8, 3, 1]
    hi = [100, 8, 5, 13, 9, 3]
    cols = [BAD, ED, PANEL, "#8a6000", PANEL, GOOD]
    y = range(len(terms))
    for i, (a, b, c) in enumerate(zip(lo, hi, cols)):
        ax.barh(i, b - a, left=a, color=c, alpha=0.75, height=0.55)
        ax.text(b + 1.5, i, "%d–%d %%" % (a, b), va="center", fontsize=9, color=INK)
    ax.set_yticks(list(y)); ax.set_yticklabels(terms, fontsize=9)
    ax.set_xlabel("typical contribution to $R_{rs}$  (%)", fontsize=10)
    ax.set_xlim(0, 118); ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.25)
    ax.set_title("Indicative magnitudes. ρ dominates everything else combined.",
                 fontsize=10.5, loc="left")

    ax2 = overlay(fig)
    box(ax2, 0.045, 0.305, 0.44, 0.215,
        "ρ is the whole game", BAD, "#fdecea",
        "Its wind dependence comes from wave facets\n"
        "sampling a NON-UNIFORM sky. That single fact\n"
        "explains the whole behaviour:\n\n"
        "  • clear sky + wind  → ρ swings widely, and no\n"
        "    bundled table covers it (Mobley 2015 is not\n"
        "    redistributable). Above ~5 m/s the software\n"
        "    REFUSES to invent a value.\n"
        "  • uniform overcast → little to sample between,\n"
        "    so ρ is stable and barely wind-dependent.\n"
        "    Overcast is the FRIENDLY case.\n\n"
        "Field rule: WHITECAPS APPEARING marks the edge\n"
        "of where ρ = 0.028 holds.", bsize=8.3)

    box(ax2, 0.515, 0.305, 0.44, 0.215,
        "What does NOT dominate", GOOD, "#eaf6ea",
        "Replicate scatter on calm water is 1–3 %. If you\n"
        "are averaging ten scans to chase that, while ρ\n"
        "carries 20–100 %, you are optimising the wrong\n"
        "term.\n\n"
        "The useful things to spend effort on, in order:\n\n"
        "  1. measure and record the WIND\n"
        "  2. get the sky-scan BEARING right (same as\n"
        "     the water scan, the commonest field error)\n"
        "  3. cross-calibrate the two instruments once\n"
        "  4. keep the panel clean, dry and LEVEL\n"
        "  5. only then, more replicates", bsize=8.3)

    box(ax2, 0.045, 0.055, 0.91, 0.225,
        "WHAT NONE OF THIS COVERS", "#8a6000", "#fff6d5",
        "These are structural, not statistical, and no amount of averaging touches "
        "them:\n\n"
        "  •  $R_{rs}$ is produced AT THE MEASUREMENT GEOMETRY, not normalised to nadir. "
        "Satellite products are\n"
        "     exact-normalised. The gap is 8–13 % at typical field geometries. Correct "
        "with the Morel f/Q ratio\n"
        "     (giop-workbench), and note that correction itself omits the air–water "
        "transmittance term.\n"
        "  •  Under OVERCAST the BRDF correction is not applicable at all: the Morel "
        "tables describe a clear-sky\n"
        "     field with a direct beam.\n"
        "  •  The panel is assumed LAMBERTIAN. Real Spectralon is close but not exact, "
        "and worse when dirty.\n"
        "  •  Everything assumes the water is HORIZONTALLY UNIFORM over the footprint "
        "and the surface is not\n"
        "     breaking. Whitecaps in the field of view invalidate the scan outright.",
        tsize=11, bsize=8.4)
    return fig


# ---------------------------------------------------------------- page 4
def page4():
    fig = plt.figure(figsize=PAGE)
    header(fig, 4, "4 · CROSS-CHECKING, AND WHAT TO DO IN THE FIELD",
           "Every check available, and what each one can and cannot catch")
    ax = overlay(fig)

    checks = [
        ("1 · $E_d$ STABILITY   (needs the irradiance sensor)", ED, "#f2edf9",
         "Bracket the station with two $E_d$ scans, before and after the target sequence.\n"
         "    ed_stability(before, after, wl)  →  worst % change + verdict\n"
         "CATCHES: the light changing mid-station, which violates the core assumption\n"
         "that $E_d$ is the same for reference and target. A 22 % drop reads 'every $R_{rs}$\n"
         "here is suspect'. With a panel alone you cannot test this without spending\n"
         "another scan, and even then the interval between is unsampled."),
        ("2 · INTER-INSTRUMENT  C(λ)   (needs panel + irradiance sensor)", GOOD, "#eaf6ea",
         "Simultaneous panel scan and $E_d$ scan → equation (F4).\n"
         "    cross_calibration_factor(wl, L_panel, Ed, ed_wl)\n"
         "CATCHES: a relative calibration offset between the two instruments, which is\n"
         "otherwise a silent multiplicative bias on every $R_{rs}$. Distinguishes a FLAT\n"
         "disagreement (gain) from a STRUCTURED one (stale calibration, tilt, shading)."),
        ("3 · ATMOSPHERIC TRANSMITTANCE   (needs the irradiance sensor)", "#8a6000",
         "#fff6d5",
         "    T(λ) = $E_d$(λ) / [ $F_0$(λ) · cos θ$_s$ ]                              (F5)\n"
         "CATCHES: an objective record of the sky state, replacing 'it looked clear'.\n"
         "Clear sky 0.6–0.8 in the visible; heavy overcast 0.1–0.3. The SHAPE matters\n"
         "too: clear-sky T rises toward the red as Rayleigh weakens, while cloud is\n"
         "nearly flat, so a FLAT T is itself evidence of cloud. In giop-workbench,\n"
         "which carries the $F_0$ table."),
        ("4 · REPLICATES", "#2c6f9b", "#eaf4f8",
         "Five water and five sky scans per station; the GUI averages and plots them\n"
         "with the individuals.\n"
         "CATCHES: random surface variability only. It does NOT catch ρ, calibration,\n"
         "or geometry error — those are systematic and shared by every replicate, so\n"
         "the scatter of replicates UNDERSTATES the real uncertainty."),
        ("5 · SANITY, on the spot", BAD, "#fdecea",
         "$R_{rs}$ positive across 400–700? Negative means over-subtraction.\n"
         "Peak where the water looks? Blue for clear, green for turbid.\n"
         "$R_{rs}$(443) of order 1e-3 to 1e-2?\n"
         "CATCHES: gross errors, immediately, while you can still re-measure. This is\n"
         "the single highest-value habit: run the GUI AT the station, not at home."),
    ]
    y = 0.895
    for title, col, fc, body in checks:
        h = 0.128
        box(ax, 0.045, y - h, 0.91, h, title, col, fc, body, tsize=10.5, bsize=8.4)
        y -= h + 0.017

    ax.text(0.5, 0.135,
            "No check here catches ρ. Nothing you can do in the field measures it. "
            "That is why the wind speed\nand the sky state must be RECORDED: they are "
            "the only route to a better ρ later.",
            fontsize=10, color=BAD, ha="center", weight="bold", va="top")
    return fig


# ---------------------------------------------------------------- page 5
def page5():
    fig = plt.figure(figsize=PAGE)
    header(fig, 5, "5 · WHAT THE ABSOLUTE SCALE BUYS",
           "Products a reflectance-only instrument cannot produce")
    ax = overlay(fig)

    # --- PAR
    ax.add_patch(FancyBboxPatch((0.045, 0.700), 0.91, 0.205,
                                boxstyle="round,pad=0.008", fc="#eaf6ea", ec=GOOD,
                                lw=2.0, transform=ax.transAxes, zorder=2))
    ax.text(0.5, 0.875, "PAR \u2014 Photosynthetically Available Radiation", fontsize=12,
            weight="bold", color=GOOD, ha="center", zorder=3)
    f = Flow(ax, 0.075, 0.845)
    f.text("A PHOTON flux, not an energy flux, because photosynthesis counts photons "
           "rather than joules:")
    f.eq(r"$\mathrm{PAR}\;=\;\int_{400}^{700} E_d(\lambda)\,"
         r"\frac{\lambda}{119.6}\;d\lambda$"
         r"$\qquad$ [$\mu$mol photons m$^{-2}$ s$^{-1}$]$\qquad$ (F6)")
    f.text("with $E_d$ in W m\u207b\u00b2 nm\u207b\u00b9 and \u03bb in nm. The constant is "
           "$N_A h c$ = 0.1196 J\u00b7m/mol, expressed\nfor nanometres and micromoles. "
           "Nothing here is fitted or tabulated \u2014 it is the photon energy "
           "$hc/\\lambda$.")
    f.gap(0.6)
    f.text("Full midday sun \u2248 2000 \u00b5mol m\u207b\u00b2 s\u207b\u00b9. Heavy overcast is one to "
           "two orders lower. Under cloud\nthat number IS the record of how much light "
           "the water column actually received, which is the\ncontext a cloudy-day "
           "spectrum needs and which no reflectance ratio can supply.")

    # --- nLw
    ax.add_patch(FancyBboxPatch((0.045, 0.475), 0.91, 0.195,
                                boxstyle="round,pad=0.008", fc="#eaf4f8", ec="#2c6f9b",
                                lw=2.0, transform=ax.transAxes, zorder=2))
    ax.text(0.5, 0.640, "nLw \u2014 normalised water-leaving radiance", fontsize=12,
            weight="bold", color="#2c6f9b", ha="center", zorder=3)
    f = Flow(ax, 0.075, 0.610)
    f.text("The quantity satellite ocean-colour products are actually distributed in:")
    f.eq(r"$nL_w(\lambda)\;=\;R_{rs}(\lambda)\;F_0(\lambda)$"
         r"$\qquad\qquad$ (F7)")
    f.text("$F_0$ is the extraterrestrial solar irradiance. This is a pure UNIT CHANGE "
           "from $R_{rs}$, carrying no\nextra assumption beyond the $F_0$ spectrum "
           "itself, and it is what makes a field spectrum\ndirectly comparable with a "
           "satellite product.")
    f.gap(0.6)
    f.text("It does NOT apply the BRDF correction, so the result is still at the "
           "measurement geometry, and\nit applies no Earth\u2013Sun distance correction "
           "($F_0$ is the mean-distance value, \u00b13.4 % annually).\n"
           "In giop-workbench, which carries the $F_0$ table.")

    # --- footprint
    ax.add_patch(FancyBboxPatch((0.045, 0.170), 0.91, 0.275,
                                boxstyle="round,pad=0.008", fc="#fff6d5", ec="#8a6000",
                                lw=2.0, transform=ax.transAxes, zorder=2))
    ax.text(0.5, 0.415, "FOOTPRINT \u2014 what patch of water a scan actually averaged over",
            fontsize=12, weight="bold", color="#8a6000", ha="center", zorder=3)
    f = Flow(ax, 0.075, 0.385)
    f.text("From the logged rangefinder distance R and the foreoptic full-angle FOV:")
    f.eq(r"across-view $=\;2R\tan(\mathrm{FOV}/2)$"
         r"$\qquad$ along-view $=\;$across$/\cos\theta_v$"
         r"$\qquad$ (F8)", size=12)
    f.text("          sensor height  =  R cos \u03b8$_v$                    "
           "spot offset from nadir  =  R sin \u03b8$_v$", size=10)
    f.gap(0.7)
    f.text("YOUR instrument, R = 3.837 m with an 8\u00b0 lens at 40\u00b0 from nadir: a "
           "0.54 \u00d7 0.70 m ellipse, 0.30 m\u00b2,\nsensor 2.94 m up, spot 2.47 m out.")
    f.gap(0.6)
    f.text("Three uses: it is the scale wave variability averages over; the 2.47 m "
           "offset says whether the\nfootprint clears YOUR OWN SHADOW; and it sets the "
           "honest comparison to satellites \u2014 an OLCI\n300 m pixel is 3.0\u00d710\u2075 "
           "times that area, a MODIS 1 km pixel 3.4\u00d710\u2076 times.")

    ax.text(0.5, 0.115,
            "All of these need the IRRADIANCE channel, not just a reflectance ratio.\n"
            "Take an $E_d$ scan at every station, even when the sky is clear.",
            fontsize=11, ha="center", color=INK, weight="bold", va="top")
    return fig


def main():
    pdf = os.path.join(OUT, "THEORY.pdf")
    figs = [page1(), page2(), page3(), page4(), page5()]
    with PdfPages(pdf) as pp:
        for f in figs:
            pp.savefig(f)
    for i, f in enumerate(figs, 1):
        f.savefig(os.path.join(OUT, "theory_p%d.png" % i), dpi=110)
        plt.close(f)
    print("wrote %s  (5 pages)" % pdf)
    print("wrote theory_p1..p5.png")


if __name__ == "__main__":
    main()
