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

    # --- SOLAR ELEVATION vs SOLAR ZENITH, as a standalone inset.
    # The instrument logs "Solar Angle", which is ELEVATION, and the two are trivially
    # confusable: they sum to 90 and both get called "the solar angle" in speech.
    # Drawn as its own small schematic in clear space rather than annotated onto the
    # main scene, where it collided with the panel.
    ix, iy, R = -2.05, 2.55, 1.05                 # inset origin and radius
    ax.add_patch(FancyBboxPatch((ix - 0.62, iy - 1.28), 2.90, 2.72,
                                boxstyle="round,pad=0.06", fc="#fffdf2",
                                ec="#d0b000", lw=1.2, zorder=2))
    ax.plot([ix, ix + 1.55], [iy, iy], color="#999", lw=1.1, zorder=3)          # horizon
    ax.plot([ix, ix], [iy, iy + 1.22], color="#999", lw=1.1, ls=":", zorder=3)  # vertical
    el = 52.0
    ex_, ey_ = ix + R * math.cos(math.radians(el)), iy + R * math.sin(math.radians(el))
    ax.annotate("", xy=(ex_, ey_), xytext=(ix, iy),
                arrowprops=dict(arrowstyle="-|>", lw=2.0, color=SUN), zorder=4)
    ax.add_patch(Circle((ex_ + 0.10, ey_ + 0.09), 0.13, fc=SUN, ec="#b8860b",
                        lw=1.1, zorder=5))
    ax.add_patch(Arc((ix, iy), 0.86, 0.86, angle=0, theta1=0, theta2=el,
                     color="#b8860b", lw=1.9, zorder=4))
    ax.text(ix + 0.52, iy + 0.13, "ELEV", fontsize=7.6, color="#8a6000",
            weight="bold", zorder=5)
    ax.add_patch(Arc((ix, iy), 1.55, 1.55, angle=0, theta1=el, theta2=90,
                     color="#666", lw=1.5, ls="--", zorder=4))
    ax.text(ix + 0.06, iy + 0.88, "ZENITH", fontsize=7.6, color="#555", zorder=5)
    ax.text(ix - 0.50, iy - 0.22,
            "ELEVATION + ZENITH = 90°\n"
            "The instrument's 'Solar Angle'\n"
            "is ELEVATION, not zenith.\n"
            "Use GPS Time, not the clock.",
            fontsize=7.8, color=INK, ha="left", va="top", zorder=5, linespacing=1.35)


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


def sky_strip(ax):
    """The first decision at any station: what sky is it? Everything else follows."""
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1)
    cases = [
        ("CLEAR / MOSTLY CLEAR", GOOD, "#eaf6ea",
         "Standard method. 40° / 135° from the sun.\n"
         "rho = 0.028 only below ~5 m/s (see WIND)."),
        ("UNIFORM OVERCAST  → GO", "#2c6f9b", "#e8f1f8",
         "NO SUN GLINT — better than clear in that one way.\n"
         "rho stable, barely wind-dependent. 135° now meaningless:\n"
         "keep 40°, pick bearing to clear the boat. NO BRDF. See page 2."),
        ("BROKEN / PATCHY CLOUD  → STOP", BAD, "#fdecea",
         "THE WORST CASE, worse than either extreme. E_d changes\n"
         "between panel and target scan; the error is unbounded and\n"
         "invisible. Wait for it to settle, or use the E_d sensor (page 2)."),
        # Quoted as 30-60 on purpose: solar ELEVATION 30-60 and solar ZENITH 30-60 are
        # the SAME window (symmetric about 45), so the elevation/zenith confusion the
        # inset warns about cannot give a wrong answer here.
        ("SUN 30–60° UP  (elev OR zenith)", "#8a6000", "#fff6d5",
         "TOO HIGH → glint, worst near 66°. NOON IS THE WORST TIME.\n"
         "TOO LOW  → weak, diffuse light. 30° = half the irradiance.\n"
         "Usable 20–70°. Best MID-MORNING and MID-AFTERNOON."),
    ]
    w = 2.47
    for i, (title, col, fc, body) in enumerate(cases):
        x = 0.10 + i * w
        ax.add_patch(FancyBboxPatch((x, 0.04), w - 0.14, 0.92,
                                    boxstyle="round,pad=0.04", fc=fc, ec=col, lw=2.0))
        ax.text(x + 0.10, 0.80, title, fontsize=9.0, weight="bold", color=col)
        ax.text(x + 0.10, 0.62, body, fontsize=6.9, color=INK, va="top",
                linespacing=1.32)
    ax.text(5.0, 1.02,
            "FIRST DECISIONS AT EVERY STATION:  WHAT SKY IS IT, AND WHERE IS THE SUN?",
            fontsize=11, weight="bold", color=INK, ha="center")


def steps(ax):
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0.18, 3.18)      # cropped so the boxes hug their text

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
        ("RECORD  ·  WIND", "#8a7000", "#fff6d5",
         "WIND sets ρ, the largest error in\n"
         "the whole measurement. NO\n"
         "ANEMOMETER NEEDED — read it\n"
         "off the water (Beaufort):\n\n"
         "  0–2   < 3 m/s    mirror, then\n"
         "                   wavelets, none\n"
         "                   breaking     ρ ok\n"
         "  3     3.5–5      crests break,\n"
         "                   SCATTERED\n"
         "                   whitecaps    limit\n"
         "  4+    > 5.5      FREQUENT\n"
         "                   whitecaps  INVALID\n\n"
         "→ WHITECAPS APPEARING = the\n"
         "  edge of where ρ = 0.028 holds.\n\n"
         "Also record: wind DIRECTION,\n"
         "a photo, cloud / sun obscured,\n"
         "panel reflectance (0.99?), and\n"
         "clear or turbid water\n"
         "  turbid → NOT nir_zero"),
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
        ax.add_patch(FancyBboxPatch((x, 0.26), w - 0.16, 2.82,
                                    boxstyle="round,pad=0.05", fc=fc, ec=col, lw=2.4))
        ax.text(x + (w - 0.16) / 2.0, 2.86, title, fontsize=12.5, weight="bold",
                color=col, ha="center")
        ax.text(x + 0.11, 2.60, body, fontsize=7.5, color=INK, va="top", ha="left",
                linespacing=1.32)


def page_two(fig):
    """Cloudy-weather working, and the products an ABSOLUTE radiometer supports."""
    gs = fig.add_gridspec(2, 2, hspace=0.13, wspace=0.09,
                          left=0.035, right=0.975, top=0.90, bottom=0.04)

    panels = [
        (gs[0, 0], "UNIFORM OVERCAST — what changes", "#2c6f9b", "#e8f1f8",
         "BETTER than clear sky:\n"
         "  • NO SUN GLINT. The specular sun beam is the dominant\n"
         "    error in above-water work. Under thick cloud there is\n"
         "    no beam, so the hardest thing to avoid goes away.\n"
         "  • rho is stable and barely wind-dependent. Its wind\n"
         "    sensitivity comes from wave facets sampling DIFFERENT\n"
         "    PARTS OF A NON-UNIFORM SKY. Under a uniform sky there\n"
         "    is little to sample between, so facet tilt stops\n"
         "    mattering. rho = 0.028 is defensible even in wind that\n"
         "    would rule it out on a clear day.\n\n"
         "WORSE:\n"
         "  • E_d much lower → raise integration time, more replicates\n"
         "  • No satellite match-up: it cannot see through cloud either\n"
         "  • DO NOT apply BRDF. The Morel f/Q tables describe a\n"
         "    clear-sky field with a direct beam; overcast is all\n"
         "    diffuse and those tables do not describe it.\n\n"
         "PROTOCOL: 135° azimuth becomes meaningless — no sun to point\n"
         "away from. KEEP 40° from nadir. Choose the bearing purely to\n"
         "clear the boat, its shadow and its wake. Record 'overcast':\n"
         "nothing in the file records the sky state for you."),

        (gs[0, 1], "BROKEN CLOUD — the worst case", BAD, "#fdecea",
         "Worse than either clear OR fully overcast.\n\n"
         "The panel method assumes E_d is THE SAME at the moment of\n"
         "the panel scan and the moment of the target scan. Under\n"
         "moving cloud it is not, and the error is unbounded and\n"
         "invisible in the output — nothing in the spectrum tells you\n"
         "the light changed between scans.\n\n"
         "OPTIONS, in order:\n"
         "  1. Wait for the sky to settle, either uniformly clear or\n"
         "     uniformly overcast.\n"
         "  2. Use the cosine-collector E_d channel (panel opposite).\n"
         "     Measuring E_d WITH the target removes the time lag.\n"
         "  3. If you must proceed: shorten the panel-to-target gap to\n"
         "     seconds, bracket each target between two panel scans,\n"
         "     take many replicates, and RECORD that the sky was\n"
         "     patchy so the data can be down-weighted later.\n\n"
         "A photograph of the sky at each station settles arguments\n"
         "later that no number can."),

        (gs[1, 0], "THE E_d SENSOR PATH — use it when light moves", "#8a6000", "#fff6d5",
         "Fit the COSINE DIFFUSER and take an irradiance scan. DARWin\n"
         "writes an 'Irr. (Target)' or 'Irr. (Ref.)' column, and the\n"
         "software uses MEASURED E_d instead of inferring it:\n\n"
         "    rrs_from_sed(water, sky, source='irradiance')\n\n"
         "Three error sources drop out at once:\n"
         "  • the panel-to-target TIME LAG  → the broken-cloud fix\n"
         "  • the PANEL REFLECTANCE entirely (no 0.99 assumption,\n"
         "    no certificate, no multiplicative bias)\n"
         "  • PANEL LEVELNESS — the collector defines the horizontal\n"
         "    plane itself\n\n"
         "Conditions: the irradiance channel must be calibrated in\n"
         "W m⁻² nm⁻¹ consistently with the radiance channel, and the\n"
         "COLLECTOR MUST BE LEVEL for the same reason the panel was.\n"
         "Its cosine error behaves better under diffuse light than\n"
         "under a low direct sun, so overcast is where it is strongest.\n\n"
         "Verified: recovers a known R_rs to 9 decimals and is exactly\n"
         "invariant to panel reflectance."),

        (gs[1, 1], "YOU HAVE A SPECTRORADIOMETER — use the absolute scale",
         GOOD, "#eaf6ea",
         "Absolute calibration supports products a reflectance-only\n"
         "instrument cannot give you.\n\n"
         "PAR — Photosynthetically Available Radiation, from measured\n"
         "E_d. A PHOTON flux, not an energy flux, because photosynthesis\n"
         "counts photons:\n\n"
         "    PAR = ∫(400–700) E_d(λ) · λ / 119.6  dλ\n"
         "    E_d in W m⁻² nm⁻¹, λ in nm → µmol m⁻² s⁻¹\n\n"
         "  Full midday sun ≈ 2000 µmol m⁻² s⁻¹. Heavy overcast is one\n"
         "  to two orders lower. Under cloud that number IS the record\n"
         "  of how much light the water column actually received.\n\n"
         "nLw — normalised water-leaving radiance, nLw = R_rs × F₀.\n"
         "  This is the quantity satellite ocean-colour products are\n"
         "  distributed in, so it is what makes your field spectrum\n"
         "  directly comparable with them.\n\n"
         "Both need the IRRADIANCE channel, not just reflectance.\n"
         "Take an E_d scan at every station even when the sky is clear."),
    ]

    for spec, title, col, fc, body in panels:
        ax = fig.add_subplot(spec)
        ax.axis("off")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.add_patch(FancyBboxPatch((0.005, 0.005), 0.99, 0.99,
                                    boxstyle="round,pad=0.008", fc=fc, ec=col, lw=2.4,
                                    transform=ax.transAxes))
        ax.text(0.5, 0.955, title, fontsize=12.5, weight="bold", color=col,
                ha="center", transform=ax.transAxes)
        ax.text(0.03, 0.90, body, fontsize=8.3, color=INK, va="top", ha="left",
                linespacing=1.34, transform=ax.transAxes)

    fig.suptitle("ABOVE-WATER R$_{rs}$  ·  page 2 — CLOUD, THE E$_d$ SENSOR, "
                 "AND ABSOLUTE PRODUCTS", fontsize=18, weight="bold", color=INK, y=0.963)
    fig.text(0.5, 0.925,
             "Overcast is a usable condition. Broken cloud is not. The irradiance "
             "channel is what makes the difference.",
             fontsize=11, color="#444", ha="center")


def main():
    fig = plt.figure(figsize=(16.5, 11.7))          # A3 landscape
    gs = fig.add_gridspec(3, 2, height_ratios=[1.42, 0.32, 1.26],
                          width_ratios=[1.5, 1.0], hspace=0.14, wspace=0.10,
                          left=0.03, right=0.985, top=0.905, bottom=0.03)
    side_view(fig.add_subplot(gs[0, 0]))
    plan_view(fig.add_subplot(gs[0, 1]))
    sky_strip(fig.add_subplot(gs[1, :]))
    steps(fig.add_subplot(gs[2, :]))

    fig.suptitle("ABOVE-WATER R$_{rs}$  —  three-scan field method   ·   "
                 "Spectral Evolution NaturaSpec Plus",
                 fontsize=19, weight="bold", color=INK, y=0.965)
    fig.text(0.5, 0.928,
             "Geometry after Mobley (1999): view 40° from nadir, 135° in azimuth from "
             "the sun.   Sky scan is the MIRROR of the water scan.",
             fontsize=11, color="#444", ha="center")

    fig2 = plt.figure(figsize=(16.5, 11.7))
    page_two(fig2)

    png = os.path.join(OUT, "FIELD_CARD.png")
    png2 = os.path.join(OUT, "FIELD_CARD_p2.png")
    pdf = os.path.join(OUT, "FIELD_CARD.pdf")
    fig.savefig(png, dpi=170)
    fig2.savefig(png2, dpi=170)
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(pdf) as pp:
        pp.savefig(fig)
        pp.savefig(fig2)
    print("wrote %s\nwrote %s\nwrote %s (2 pages — print double-sided)"
          % (png, png2, pdf))


if __name__ == "__main__":
    main()
