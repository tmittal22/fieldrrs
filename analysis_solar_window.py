"""Where does the sun have to be? Compute the window instead of asserting it.

    python analysis_solar_window.py

`solar.py` shipped SOLAR_ZENITH_MIN/MAX = 20/60 deg with the justification "either the
sun is high enough that glint is hard to avoid, or low enough that the signal is weak
and the atmospheric path long". Both halves are plausible and neither was ever computed.
This computes them.

Three independent constraints, each a function of solar zenith theta_s, at the shipped
viewing geometry (40 deg from nadir, 135 deg in azimuth from the sun):

  1 DIRECT SUN GLINT.  A wave facet reflects the sun specularly into the sensor when its
    normal bisects the sun and view directions. That required tilt beta(theta_s) is exact
    geometry. Whether such a facet EXISTS is Cox & Munk (1954): sea-surface slopes are
    near-Gaussian with variance sigma^2 = 0.003 + 0.00512 W. So the glint risk is set by
    tan(beta)/sigma, in units of standard deviations of slope. Small = dangerous.

  2 SIGNAL.  E_d on the horizontal falls as cos(theta_s), and the atmospheric path grows
    as the airmass ~ 1/cos(theta_s), so the direct beam is attenuated further. Low sun is
    photon-starved AND increasingly diffuse.

  3 PLATFORM SHADOW.  The operator's shadow falls on the anti-solar bearing with length
    h*tan(theta_s), which grows without bound as the sun drops.

RESULT, and two of the three surprised me:

  * Glint is WORST AT HIGH SUN, not low. The required facet tilt is MINIMISED at about
    65 deg elevation (15.2 deg), so the shipped docstring's "high enough that glint is
    hard to avoid" is right and is now quantified. Dropping from 60 to 30 deg elevation
    cuts the Cox-Munk glint weight by 15x (0.254 -> 0.017 at 5 m/s).
  * SHADOW DOES NOT BIND. Reach is the wrong test: the shadow runs along the anti-solar
    bearing while the spot sits 45 deg off it, so the miss distance is a FIXED 1.74 m
    lateral, independent of solar elevation. A person never shadows the spot at this
    standoff. A boat or a pier is a different object and can.
  * Only SIGNAL argues for high sun, and it is the reason the window has a bottom.

So the two ends are set by glint (top) and signal (bottom), and the shadow term that
motivates IOCCG's preference for 90 deg azimuth is about the STRUCTURE, not the operator.

Cox & Munk (1954), doi:10.1364/JOSA.44.000838.
"""

import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

VIEW_ZENITH = 40.0
REL_AZIMUTH = 135.0
ALT_AZIMUTH = 90.0


def facet_tilt_deg(theta_s_deg, theta_v_deg=VIEW_ZENITH, dphi_deg=REL_AZIMUTH):
    """Facet tilt that specularly reflects the sun into the sensor, degrees.

    The facet normal must bisect the direction TO the sun and the direction TO the
    sensor, both taken as unit vectors pointing up from the surface:

        n  proportional to  s_hat + v_hat
        cos(beta) = (cos theta_s + cos theta_v) / |s_hat + v_hat|

    Exact for a flat facet; no wave statistics enter here.
    """
    ts, tv = math.radians(theta_s_deg), math.radians(theta_v_deg)
    dphi = math.radians(dphi_deg)
    sdotv = math.sin(ts) * math.sin(tv) * math.cos(dphi) + math.cos(ts) * math.cos(tv)
    norm = math.sqrt(2.0 + 2.0 * sdotv)
    return math.degrees(math.acos((math.cos(ts) + math.cos(tv)) / norm))


def cox_munk_sigma(wind_ms):
    """RMS sea-surface slope (dimensionless), Cox & Munk (1954) clean-surface fit."""
    return math.sqrt(0.003 + 0.00512 * wind_ms)


def glint_sigmas(theta_s_deg, wind_ms, **kw):
    """How many slope standard deviations away the glint-producing facet is.

    Below ~2 sigma such facets are common and direct glint is likely in the field of
    view; beyond ~3 sigma they are rare.
    """
    return math.tan(math.radians(facet_tilt_deg(theta_s_deg, **kw))) / cox_munk_sigma(wind_ms)


def shadow_reach_m(theta_s_deg, operator_h=2.0):
    """How far the operator's shadow extends along the anti-solar bearing."""
    return operator_h * math.tan(math.radians(theta_s_deg))


def shadow_clearance_m(range_m=3.837, theta_v_deg=VIEW_ZENITH, dphi_deg=REL_AZIMUTH):
    """Perpendicular distance from the target spot to the shadow AXIS.

    Reach alone is the wrong test. The shadow runs along the anti-solar bearing; the
    spot sits at 180 - dphi degrees away from that bearing, so the miss distance is a
    fixed lateral offset INDEPENDENT of solar elevation. A shadow that reaches 11 m
    still misses if it passes 1.75 m to the side.
    """
    off_axis = math.radians(180.0 - dphi_deg)
    return range_m * math.sin(math.radians(theta_v_deg)) * math.sin(off_axis)


def glint_weight(theta_s_deg, wind_ms, **kw):
    """Cox-Munk Gaussian weight of the glint-producing facet, relative to a flat sea.

    exp(-tan^2(beta) / (2 sigma^2)). This is the shape that makes rho depend on solar
    zenith and wind at all; it is not an extra error on top of rho, it is WHY a single
    rho = 0.028 cannot be right everywhere.
    """
    t = math.tan(math.radians(facet_tilt_deg(theta_s_deg, **kw)))
    return math.exp(-t * t / (2.0 * cox_munk_sigma(wind_ms) ** 2))


def spot_offset_m(range_m=3.837, theta_v_deg=VIEW_ZENITH):
    return range_m * math.sin(math.radians(theta_v_deg))


def main():
    zen = [z + 0.5 for z in range(0, 85)]
    elev = [90.0 - z for z in zen]

    print("theta_s  elev   beta(135)  n_sig@5  glintwt@5  glintwt@10  Ed/Ed(0)  "
          "shadow(m)")
    for z in (0, 10, 20, 25, 30, 40, 50, 55, 60, 65, 70, 75, 80):
        print("%5.0f  %5.0f   %7.1f  %7.2f  %9.3f  %10.3f  %8.2f  %7.2f"
              % (z, 90 - z, facet_tilt_deg(z), glint_sigmas(z, 5.0),
                 glint_weight(z, 5.0), glint_weight(z, 10.0),
                 math.cos(math.radians(z)), shadow_reach_m(z)))
    print("\nShadow LATERAL clearance at 135 deg azimuth, 3.837 m range: %.2f m "
          "(independent of solar elevation)." % shadow_clearance_m())
    print("At 90 deg azimuth: %.2f m." % shadow_clearance_m(dphi_deg=ALT_AZIMUTH))

    b135 = [facet_tilt_deg(z) for z in zen]
    b90 = [facet_tilt_deg(z, dphi_deg=ALT_AZIMUTH) for z in zen]
    s5 = [glint_sigmas(z, 5.0) for z in zen]
    s10 = [glint_sigmas(z, 10.0) for z in zen]
    ed = [math.cos(math.radians(z)) for z in zen]
    shad = [shadow_reach_m(z) for z in zen]
    spot = spot_offset_m()

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0))

    ax = axes[0]
    ax.plot(elev, b135, lw=2.2, color="#1f7a99", label="135° azimuth (shipped)")
    ax.plot(elev, b90, lw=1.6, color="#6a3d9a", ls="--", label="90° azimuth (alt)")
    i = min(range(len(zen)), key=lambda k: b135[k])
    ax.plot(elev[i], b135[i], "o", color="#c0392b", ms=8, zorder=5)
    ax.annotate("worst case\n%.0f° elevation\nfacet tilt %.1f°" % (elev[i], b135[i]),
                xy=(elev[i], b135[i]), xytext=(elev[i] - 34, b135[i] + 9),
                fontsize=9, color="#c0392b",
                arrowprops=dict(arrowstyle="->", color="#c0392b"))
    ax.set_xlabel("solar elevation (deg above horizon)")
    ax.set_ylabel("facet tilt needed to put the sun in the FOV (deg)")
    ax.set_title("1 · Direct glint geometry\nLOWER = more dangerous", fontsize=10.5,
                 loc="left")
    ax.legend(fontsize=8.5); ax.grid(alpha=0.3); ax.set_xlim(0, 90)

    ax = axes[1]
    ax.plot(elev, s5, lw=2.2, color="#2e7d32", label="wind 5 m/s")
    ax.plot(elev, s10, lw=2.2, color="#d9534f", label="wind 10 m/s")
    ax.axhspan(0, 2, color="#c0392b", alpha=0.13)
    ax.axhline(2, color="#c0392b", lw=1.2, ls=":")
    ax.text(3, 1.0, "glint-producing facets COMMON\n(< 2 sigma of slope)", fontsize=8.5,
            color="#c0392b", va="center")
    ax.set_xlabel("solar elevation (deg above horizon)")
    ax.set_ylabel("glint facet tilt / RMS slope  (sigma)")
    ax.set_title("2 · Does such a facet exist?  Cox & Munk (1954)\nLOWER = more "
                 "dangerous", fontsize=10.5, loc="left")
    ax.legend(fontsize=8.5); ax.grid(alpha=0.3); ax.set_xlim(0, 90); ax.set_ylim(0, 6)

    ax = axes[2]
    ax.plot(elev, ed, lw=2.2, color="#f5a623", label="$E_d$ / $E_d$(overhead)")
    ax.set_xlabel("solar elevation (deg above horizon)")
    ax.set_ylabel("relative $E_d$ on the horizontal", color="#b07800")
    ax.set_ylim(0, 1.05); ax.grid(alpha=0.3); ax.set_xlim(0, 90)
    ax2 = ax.twinx()
    ax2.plot(elev, shad, lw=2.0, color="#1a1a1a", ls="--",
             label="operator shadow REACH (2 m)")
    clr = shadow_clearance_m()
    ax2.axhline(clr, color="#1f7a99", lw=2.0)
    ax2.text(7, 0.42, "lateral CLEARANCE %.2f m — constant, never crossed" % clr,
             fontsize=8.5, color="#1f7a99", weight="bold")
    ax2.set_ylabel("distance (m)"); ax2.set_ylim(0, 12)
    ax.set_title("3 · Signal falls and the shadow grows\nas the sun drops",
                 fontsize=10.5, loc="left")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8.5, loc="center right")

    for a in axes:
        a.axvspan(30, 60, color="#2e7d32", alpha=0.10, zorder=0)
    fig.suptitle(
        "Where the sun has to be — computed, at 40° view zenith / 135° relative "
        "azimuth.   Green band = elevation 30–60° (= zenith 30–60°, symmetric about "
        "45°).\nIt is the best available COMPROMISE, not a glint-free zone: glint is "
        "never eliminated at moderate sun, which is exactly what ρ = 0.028 exists to "
        "absorb.", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig("figures_solar_window.png", dpi=130)
    print("\nwrote figures_solar_window.png")

    print("\n--- candidate windows, all quoted as SOLAR ELEVATION ---")
    for name, lo_e, hi_e in (("shipped  (zenith 20-60) = elevation", 30.0, 70.0),
                             ("proposed              elevation", 20.0, 60.0),
                             ("intersection          elevation", 30.0, 60.0)):
        print("%-38s %2.0f-%2.0f   worst glint wt %.3f @5 m/s, min Ed %.2f"
              % (name, lo_e, hi_e,
                 max(glint_weight(90 - e, 5.0) for e in (lo_e, hi_e)),
                 math.sin(math.radians(lo_e))))
    print("\nNOTE elevation 30-60 == zenith 30-60: the window is SYMMETRIC about 45 deg,")
    print("so quoting '30 to 60 degrees' is safe whichever convention the reader assumes.")


if __name__ == "__main__":
    main()
