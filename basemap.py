"""Fetch a web-map tile mosaic for a lat/lon box. No new dependencies.

Uses the standard slippy-map (Web Mercator, EPSG:3857) tile scheme, so it works with any
XYZ provider. Tiles are cached under .tilecache/ so a re-run is offline and instant.

    img, extent = mosaic(lat_min, lat_max, lon_min, lon_max, zoom=16)
    ax.imshow(img, extent=extent, origin="upper")

`extent` is returned in DEGREES for convenience, which is only valid over a small box:
Web Mercator's y is nonlinear in latitude, and pretending otherwise stretches the image.
Over ~1 km at 67 N the error is well under a pixel, and `check_distortion` measures it
rather than assuming it.

Imagery is Esri World Imagery by default. Attribution is REQUIRED when you publish a
figure using it; `ATTRIBUTION[provider]` is the string to print.
"""

import io
import math
import os

from PIL import Image

PROVIDERS = {
    "esri_imagery": "https://server.arcgisonline.com/ArcGIS/rest/services/"
                    "World_Imagery/MapServer/tile/{z}/{y}/{x}",
    "esri_topo": "https://server.arcgisonline.com/ArcGIS/rest/services/"
                 "World_Topo_Map/MapServer/tile/{z}/{y}/{x}",
    "osm": "https://tile.openstreetmap.org/{z}/{x}/{y}.png",
}

ATTRIBUTION = {
    "esri_imagery": "Imagery: Esri, Maxar, Earthstar Geographics, and the GIS User "
                    "Community",
    "esri_topo": "Esri, HERE, Garmin, USGS, NGA",
    "osm": "(c) OpenStreetMap contributors",
}

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".tilecache")


def deg2num(lat, lon, z):
    n = 2.0 ** z
    x = (lon + 180.0) / 360.0 * n
    la = math.radians(lat)
    y = (1.0 - math.asinh(math.tan(la)) / math.pi) / 2.0 * n
    return x, y


def num2deg(x, y, z):
    n = 2.0 ** z
    lon = x / n * 360.0 - 180.0
    lat = math.degrees(math.atan(math.sinh(math.pi * (1.0 - 2.0 * y / n))))
    return lat, lon


def _tile(provider, z, x, y, timeout=20):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, "%s_%d_%d_%d.png" % (provider, z, x, y))
    if os.path.exists(path):
        return Image.open(path).convert("RGB")
    import requests
    url = PROVIDERS[provider].format(z=z, x=x, y=y)
    r = requests.get(url, timeout=timeout,
                     headers={"User-Agent": "fieldrrs/1.0 (research; contact "
                                            "mittal.tushar22@gmail.com)"})
    r.raise_for_status()
    im = Image.open(io.BytesIO(r.content)).convert("RGB")
    im.save(path)
    return im


def mosaic(lat_min, lat_max, lon_min, lon_max, zoom=16, provider="esri_imagery"):
    """Return (PIL image, (lon_min, lon_max, lat_min, lat_max)) covering the box."""
    x0, y0 = deg2num(lat_max, lon_min, zoom)      # NW corner -> smallest x, smallest y
    x1, y1 = deg2num(lat_min, lon_max, zoom)
    xi0, yi0, xi1, yi1 = (int(math.floor(x0)), int(math.floor(y0)),
                          int(math.floor(x1)), int(math.floor(y1)))
    w, h = (xi1 - xi0 + 1), (yi1 - yi0 + 1)
    out = Image.new("RGB", (w * 256, h * 256))
    for i in range(w):
        for j in range(h):
            try:
                out.paste(_tile(provider, zoom, xi0 + i, yi0 + j), (i * 256, j * 256))
            except Exception:
                pass                                   # a missing tile stays black
    lat_n, lon_w = num2deg(xi0, yi0, zoom)
    lat_s, lon_e = num2deg(xi1 + 1, yi1 + 1, zoom)
    return out, (lon_w, lon_e, lat_s, lat_n)


def check_distortion(extent, npix_y):
    """Max latitude error, in pixels, from plotting Mercator tiles on a linear axis."""
    lon_w, lon_e, lat_s, lat_n = extent
    z = 20
    ytop = deg2num(lat_n, lon_w, z)[1]
    ybot = deg2num(lat_s, lon_w, z)[1]
    worst = 0.0
    for k in range(1, 20):
        f = k / 20.0
        lat_lin = lat_n + f * (lat_s - lat_n)              # linear in degrees
        lat_mer = num2deg(0, ytop + f * (ybot - ytop), z)[0]   # true Mercator position
        worst = max(worst, abs(lat_lin - lat_mer) / abs(lat_n - lat_s))
    return worst * npix_y
