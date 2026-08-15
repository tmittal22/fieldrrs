"""Generate the printable field card (PNG + PDF).

This is a DEVELOPMENT script and is the only thing in the repository that needs
matplotlib. The field package itself stays pure standard library; the card it produces
is a static image you print and take with you.

    python make_field_card.py
"""

import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, FancyArrow, FancyBboxPatch, Polygon, Wedge

OUT = os.path.dirname(os.path.abspath(__file__))

SUN = "#f5a623"
WATER = "#1f7a99"
SKY = "#7fb3d5"
PANEL = "#d9534f"
INSTR = "#333333"
GOOD = "#2e7d32"
BAD = "#c0392b"
INK = "#1a1a1a"

VIEW_ZENITH = 40.0        # from nadir
REL_AZ = 135.0


def instrument(ax, x, y, angle_deg, scale=1.0, label=None):
    """A stylised NaturaSpec Plus: body, pistol grip, foreoptic, sighting the target."""
    a = math.radians(angle_deg)
    ca, sa = math.cos(a), math.sin(a)

    def R(px, py):
        return (x + (px * ca - py * sa) * scale, y + (px * sa + py * ca) * scale)

    body = Polygon([R(-0.55, -0.16), R(0.40, -0.16), R(0.40, 0.16), R(-0.55, 0.16)],
                   closed=True, fc=INSTR, ec="black", lw=1.2, zorder=6)
    ax.add_patch(body)
    lens = Polygon([R(0.40, -0.13), R(0.66, -0.07), R(0.66, 0.07), R(0.40, 0.13)],
                   closed=True, fc="#888", ec="black", lw=1.0, zorder=6)
    ax.add_patch(lens)
    grip = Polygon([R(-0.28, -0.16), R(-0.06, -0.16), R(-0.10, -0.52), R(-0.30, -0.52)],
                   closed=True, fc="#555", ec="black", lw=1.0, zorder=5)
    ax.add_patch(grip)
    sx, sy = R(-0.20, 0.16)
    ax.add_patch(Circle((sx, sy), 0.055 * scale, fc="#0d0", ec="black", lw=0.8, zorder=7))
    if label:
        lx, ly = R(-0.75, 0.42)
        ax.text(lx, ly, label, fontsize=8.5, color=INK, ha="center", weight="bold")


def side_view(ax):
    ax.set_title("A.  SIDE VIEW  —  the three scans and their angles",
                 fontsize=13, weight="bold", loc="left", color=INK)
    ax.set_xlim(-5.2, 5.6)
    ax.set_ylim(-2.5, 4.6)
    ax.axis("off")

    # water surface
    ax.fill_between([-5.2, 5.6], -2.5, 0, color=WATER, alpha=0.30, zorder=0)
    ax.plot([-5.2, 5.6], [0, 0], color=WATER, lw=2.5, zorder=1)
    ax.text(5.45, -0.35, "water", fontsize=10, color=WATER, ha="right", weight="bold")

    ox, oy = 0.0, 0.9          # operator / instrument origin
    ax.plot([ox, ox], [0, 3.9], color="#bbb", lw=1.0, ls=":", zorder=1)
    ax.text(ox + 0.08, 3.95, "ZENITH (straight up)", fontsize=8, color="#666")
    ax.plot([ox, ox], [0, -2.3], color="#bbb", lw=1.0, ls=":", zorder=1)
    ax.text(ox + 0.08, -2.35, "NADIR (straight down)", fontsize=8, color="#666")
    ax.plot([-5.0, 5.4], [oy, oy], color="#bbb", lw=1.0, ls="--", zorder=1)
    ax.text(-4.95, oy + 0.12, "horizontal", fontsize=8, color="#666")

    # sun
    ax.add_patch(Circle((-3.7, 3.6), 0.42, fc=SUN, ec="#b8860b", lw=1.5, zorder=5))
    for k in range(12):
        th = k * math.pi / 6
        ax.plot([-3.7 + 0.52 * math.cos(th), -3.7 + 0.72 * math.cos(th)],
                [3.6 + 0.52 * math.sin(th), 3.6 + 0.72 * math.sin(th)],
                color=SUN, lw=2, zorder=5)
    ax.text(-3.7, 4.25, "SUN", fontsize=11, weight="bold", color="#b8860b", ha="center")

    # --- WATER ray: 40 deg from nadir  =  50 deg below horizontal
    ang = math.radians(90 - VIEW_ZENITH)     # below horizontal
    L = 2.9
    wx, wy = ox + L * math.cos(-ang), oy + L * math.sin(-ang)
    ax.annotate("", xy=(wx, wy), xytext=(ox, oy),
                arrowprops=dict(arrowstyle="-|>", lw=2.6, color=WATER,
                                shrinkA=0, shrinkB=0), zorder=4)
    ax.add_patch(Arc((ox, oy), 2.6, 2.6, angle=0, theta1=-90, theta2=-(90 - VIEW_ZENITH),
                     color=WATER, lw=1.6, zorder=3))
    ax.text(ox + 0.62, oy - 1.55, "%.0f°\nfrom nadir" % VIEW_ZENITH, fontsize=9.5,
            color=WATER, ha="center", weight="bold")
    ax.add_patch(Arc((ox, oy), 3.9, 3.9, angle=0, theta1=-(90 - VIEW_ZENITH), theta2=0,
                     color="#999", lw=1.2, ls="--", zorder=3))
    ax.text(ox + 2.15, oy - 0.62, "= %.0f° below\n   horizontal" % (90 - VIEW_ZENITH),
            fontsize=9, color="#555", ha="left")
    ax.text(wx + 0.15, wy - 0.28, "3.  WATER  L_t", fontsize=11, weight="bold",
            color=WATER)

    # --- SKY ray: 40 deg from zenith = 50 deg above horizontal, SAME bearing
    sx2, sy2 = ox + L * math.cos(ang), oy + L * math.sin(ang)
    ax.annotate("", xy=(sx2, sy2), xytext=(ox, oy),
                arrowprops=dict(arrowstyle="-|>", lw=2.6, color=SKY,
                                shrinkA=0, shrinkB=0), zorder=4)
    ax.add_patch(Arc((ox, oy), 2.6, 2.6, angle=0, theta1=90 - VIEW_ZENITH, theta2=90,
                     color=SKY, lw=1.6, zorder=3))
    ax.text(ox + 0.60, oy + 1.55, "%.0f°\nfrom zenith" % VIEW_ZENITH, fontsize=9.5,
            color="#2c6f9b", ha="center", weight="bold")
    ax.text(sx2 + 0.15, sy2 + 0.18, "2.  SKY  L_sky", fontsize=11, weight="bold",
            color="#2c6f9b")

    # the mirror relationship
    ax.text(ox + 3.15, oy + 0.12,
            "SKY and WATER are MIRRORS\nthrough the horizontal:\nSAME compass bearing",
            fontsize=8.6, color=INK, ha="left",
            bbox=dict(boxstyle="round,pad=0.35", fc="#fff6d5", ec="#d0b000", lw=1))

    instrument(ax, ox, oy, -(90 - VIEW_ZENITH), 0.95)

    # --- PANEL
    px, py = -3.5, 0.45
    ax.add_patch(FancyBboxPatch((px - 0.62, py - 0.09), 1.24, 0.18,
                                boxstyle="round,pad=0.02", fc="white", ec=PANEL,
                                lw=2.2, zorder=6))
    ax.plot([px - 0.62, px + 0.62], [py - 0.11, py - 0.11], color="#999", lw=1)
    ax.plot([px, px], [py - 0.5, py - 0.11], color="#999", lw=2)
    ax.text(px, py + 0.30, "1.  PANEL  L_p", fontsize=11, weight="bold", color=PANEL,
            ha="center")
    ax.text(px, py - 0.85, "LEVEL, face up\nview from straight down",
            fontsize=8.6, color=PANEL, ha="center")
    ax.annotate("", xy=(px, py + 0.16), xytext=(px, py + 1.35),
                arrowprops=dict(arrowstyle="-|>", lw=2.2, color=PANEL), zorder=4)
    instrument(ax, px, py + 1.75, -90, 0.75)

    ax.text(-5.0, -1.55,
            "E$_d$ = π · L$_p$ / R$_{panel}$        "
            "L$_w$ = L$_t$ − ρ · L$_{sky}$        "
            "R$_{rs}$ = L$_w$ / E$_d$",
            fontsize=12, color=INK, family="serif",
            bbox=dict(boxstyle="round,pad=0.45", fc="#eef5ff", ec="#8aa8c8", lw=1.2))


def plan_view(ax):
    ax.set_title("B.  PLAN VIEW  —  which way to point (looking down)",
                 fontsize=13, weight="bold", loc="left", color=INK)
    ax.set_xlim(-1.75, 1.75)
    ax.set_ylim(-1.75, 1.60)
    ax.set_aspect("equal")
    ax.axis("off")

    sun_az = 160.0          # worked example, chosen clear of the cardinal labels
    R = 1.0

    def pt(az, r=R):
        a = math.radians(90 - az)
        return r * math.cos(a), r * math.sin(a)

    # glint zone: within +/-45 deg of the sun bearing
    ax.add_patch(Wedge((0, 0), R, 90 - (sun_az + 45), 90 - (sun_az - 45),
                       fc=BAD, alpha=0.15, zorder=1))
    ax.text(*pt(sun_az, 0.62), s="GLINT\navoid", fontsize=8.5, color=BAD,
            ha="center", va="center", weight="bold", zorder=6,
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.75))

    ax.add_patch(Circle((0, 0), R, fill=False, ec="#999", lw=1.4, zorder=2))
    # cardinals INSIDE the ring, so nothing outside can collide with them
    for az, lab in ((0, "N"), (90, "E"), (180, "S"), (270, "W")):
        ax.text(*pt(az, 0.87), s=lab, fontsize=11, weight="bold", ha="center",
                va="center", color="#777", zorder=3)
    for az in range(0, 360, 30):
        ax.plot(*zip(pt(az, 0.97), pt(az, 1.0)), color="#bbb", lw=1)

    # sun direction
    ax.annotate("", xy=pt(sun_az), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", lw=3, color=SUN), zorder=4)
    ax.add_patch(Circle(pt(sun_az, 1.12), 0.11, fc=SUN, ec="#b8860b", lw=1.3, zorder=5))
    ax.text(*pt(sun_az, 1.46), s="SUN  %.0f°" % sun_az, fontsize=10.5, weight="bold",
            color="#b8860b", ha="center", va="center")

    # the two acceptable view bearings
    for sign in (-1, +1):
        az = (sun_az + sign * REL_AZ) % 360
        ax.annotate("", xy=pt(az), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", lw=3, color=GOOD), zorder=4)
        ax.text(*pt(az, 1.34), s="POINT HERE\n%.0f°" % az, fontsize=10, weight="bold",
                color=GOOD, ha="center", va="center")
        a0, a1 = 90 - sun_az, 90 - az
        ax.add_patch(Arc((0, 0), 1.30, 1.30, angle=0,
                         theta1=min(a0, a1), theta2=max(a0, a1),
                         color=GOOD, lw=1.3, ls="--", zorder=3))
        ax.text(*pt(sun_az + sign * REL_AZ / 2.0, 0.72), s="%.0f°" % REL_AZ,
                fontsize=9, color=GOOD, ha="center", va="center", weight="bold",
                bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none", alpha=0.8))

    ax.add_patch(Circle((0, 0), 0.05, fc=INK, zorder=6))
    ax.text(0, -1.62,
            "Both bearings are %.0f° in azimuth from the sun. Pick whichever is clear\n"
            "of the boat, its shadow and its wake.   The GUI computes these for you\n"
            "from latitude, longitude, date and time." % REL_AZ,
            fontsize=9.2, color=INK, ha="center", va="center")


def steps(ax):
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0.62, 3.18)      # cropped so the boxes hug their text

    boxes = [
        ("1.  PANEL", PANEL, "white",
         "White reference panel (Spectralon),\n"
         "LEVEL and face up, in full sun.\n"
         "No shadow on it, yours included.\n"
         "View from straight down.\n"
         "Do not touch the white face.\n\n"
         "DARWin: take as the REFERENCE.\n"
         "Re-take every ~10 min, and the\n"
         "moment the light changes."),
        ("2.  SKY", "#2c6f9b", "white",
         "Sky at 40° from ZENITH\n"
         "= 50° ABOVE horizontal.\n\n"
         "SAME compass bearing as the\n"
         "water scan. Not the opposite\n"
         "one. This is the single most\n"
         "common field mistake.\n\n"
         "DARWin: TARGET scan. Save."),
        ("3.  WATER", WATER, "white",
         "Water at 40° from NADIR\n"
         "= 50° BELOW horizontal,\n"
         "135° in azimuth from the sun.\n"
         "FROM A BOAT OR PIER use 90°:\n"
         "135° looks back at the hull\n"
         "or its shadow (IOCCG v3.0).\n"
         "No glint, foam or wake.\n\n"
         "DARWin: TARGET scan. Save.\n"
         "Repeat 2 and 3 five times."),
        ("RECORD", "#8a7000", "#fff6d5",
         "WIND SPEED (m/s) — sets ρ, the\n"
         "  largest error in the whole\n"
         "  measurement. ρ = 0.028 only\n"
         "  valid below ~5 m/s.\n"
         "Wind DIRECTION + a photo\n"
         "Cloud, or sun obscured?\n"
         "Water clear or turbid?\n"
         "  turbid → do NOT use nir_zero\n"
         "Panel reflectance used (0.99?)\n\n"
         "The instrument already logs\n"
         "Range, Tilt, GPS Time and\n"
         "Solar Angle (= ELEVATION,\n"
         "not zenith — verified)."),
        ("CHECK ON THE SPOT", GOOD, "#eaf6ea",
         "⚠ THE CLOCK MAY BE WRONG.\n"
         "  Use GPS Time for solar\n"
         "  geometry. On a real file the\n"
         "  clock read 05:56 with the\n"
         "  sun below the horizon.\n\n"
         "R_rs positive across 400-700?\n"
         "  negative = over-subtracted\n"
         "Peak where the water looks?\n"
         "Replicates on top of each other?\n"
         "R_rs(443) order 1e-3 to 1e-2?\n\n"
         "FOV × range = footprint. An 8°\n"
         "lens at 3.8 m sees ~0.5 m of\n"
         "water; the GUI plots it.\n\n"
         "Copy the .sed files off tonight."),
    ]
    n = len(boxes)
    w = 9.76 / n
    for i, (title, col, fc) in enumerate([(b[0], b[1], b[2]) for b in boxes]):
        body = boxes[i][3]
        x = 0.12 + i * w
        ax.add_patch(FancyBboxPatch((x, 0.70), w - 0.16, 2.38,
                                    boxstyle="round,pad=0.05", fc=fc, ec=col, lw=2.4))
        ax.text(x + (w - 0.16) / 2.0, 2.86, title, fontsize=12.5, weight="bold",
                color=col, ha="center")
        ax.text(x + 0.11, 2.58, body, fontsize=8.6, color=INK, va="top", ha="left",
                linespacing=1.32)


def main():
    fig = plt.figure(figsize=(16.5, 11.7))          # A3 landscape
    gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 0.92], width_ratios=[1.5, 1.0],
                          hspace=0.16, wspace=0.10,
                          left=0.03, right=0.985, top=0.905, bottom=0.03)
    side_view(fig.add_subplot(gs[0, 0]))
    plan_view(fig.add_subplot(gs[0, 1]))
    steps(fig.add_subplot(gs[1, :]))

    fig.suptitle("ABOVE-WATER R$_{rs}$  —  three-scan field method   ·   "
                 "Spectral Evolution NaturaSpec Plus",
                 fontsize=19, weight="bold", color=INK, y=0.965)
    fig.text(0.5, 0.928,
             "Geometry after Mobley (1999): view 40° from nadir, 135° in azimuth from "
             "the sun.   Sky scan is the MIRROR of the water scan.",
             fontsize=11, color="#444", ha="center")

    png = os.path.join(OUT, "FIELD_CARD.png")
    pdf = os.path.join(OUT, "FIELD_CARD.pdf")
    fig.savefig(png, dpi=170)
    fig.savefig(pdf)
    print("wrote %s\nwrote %s" % (png, pdf))


if __name__ == "__main__":
    main()
