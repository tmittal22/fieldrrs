"""Spectral Evolution ``.sed`` reader. Pure standard library, no numpy.

Column vocabulary taken from the string table of ``DARWin2.exe`` in the NaturaSpec Plus
installer, not from a third-party parser: DARWin writes ``Irr. (Ref.)``, where the
parsers in circulation expect ``Irrad. (Ref.)`` and silently drop the column.

File layout::

    Version: 2.1
    Instrument: ...
    Comment: station 3 water
    Data:
    Wvl<TAB>Rad. (Ref.)<TAB>Rad. (Target)<TAB>Reflect. %
    350.0<TAB>1.234e+000<TAB>5.678e-002<TAB>4.601
"""

from __future__ import annotations

import os
import re

__all__ = ["SedSpectrum", "read_sed", "read_folder", "COLUMN_ALIASES"]

COLUMN_ALIASES = {
    "wavelength": ("Wvl",),
    "rad_ref": ("Rad. (Ref.)",),
    "rad_target": ("Rad. (Target)",),
    "irr_ref": ("Irr. (Ref.)", "Irrad. (Ref.)"),
    "irr_target": ("Irr. (Target)", "Irrad. (Target)"),
    "dn_ref": ("DN (Ref.)", "Ref. DN"),
    "dn_target": ("DN (Target)", "Tgt. DN"),
    "reflectance": ("Reflect. %", "Reflect. [1.0]", "Reflect."),
}

_DATA = re.compile(r"^\s*Data:\s*$", re.IGNORECASE)


class SedSpectrum(object):
    """One scan. ``columns`` maps canonical name -> list of floats."""

    def __init__(self, path, header, columns, raw_columns, reflectance_scale):
        self.path = path
        self.name = os.path.basename(path)
        self.header = header
        self.columns = columns
        self.raw_columns = raw_columns
        self.reflectance_scale = reflectance_scale

    # -- data ---------------------------------------------------------------
    @property
    def wavelength(self):
        return self.columns["wavelength"]

    @property
    def reflectance(self):
        """Target/reference ratio as a FRACTION, never percent."""
        if "reflectance" not in self.columns:
            raise KeyError(
                "%s has no reflectance column; found %s"
                % (self.name, sorted(self.raw_columns))
            )
        s = self.reflectance_scale or 1.0
        return [v / s for v in self.columns["reflectance"]]

    def need(self, key):
        if key not in self.columns:
            raise KeyError(
                "%s has no '%s' column (expected one of %s). Found: %s.\n"
                "Re-export from DARWin with radiance columns enabled, or use the "
                "reflectance path." % (self.name, key, COLUMN_ALIASES[key],
                                       sorted(self.raw_columns))
            )
        return self.columns[key]

    @property
    def radiance_target(self):
        return self.need("rad_target")

    @property
    def radiance_reference(self):
        return self.need("rad_ref")

    def has(self, key):
        return key in self.columns

    # -- metadata -----------------------------------------------------------
    @property
    def comment(self):
        return self.header.get("Comment", "")

    @property
    def when(self):
        return ("%s %s" % (self.header.get("Date", ""),
                           self.header.get("Time", ""))).strip()

    @property
    def latitude(self):
        return _coord(self.header.get("Latitude"))

    @property
    def longitude(self):
        return _coord(self.header.get("Longitude"))

    def clip(self, lo, hi):
        keep = [i for i, w in enumerate(self.wavelength) if lo <= w <= hi]
        cols = {k: [v[i] for i in keep] for k, v in self.columns.items()}
        raw = {k: [v[i] for i in keep] for k, v in self.raw_columns.items()}
        return SedSpectrum(self.path, self.header, cols, raw, self.reflectance_scale)

    def __repr__(self):
        w = self.wavelength
        return "SedSpectrum(%s, %d bands %.1f-%.1f nm, %s)" % (
            self.name, len(w), w[0], w[-1], sorted(self.columns))


#: Columns DARWin can write that fieldrrs does not need. Listed so the diagnostic can
#: report "known, not used" rather than "unrecognised", which would read as a defect.
KNOWN_UNUSED = ("Tgt./Ref. %", "Tgt./Ref.", "Chan#", "Counts",
                "Raw Detector Counts", "Normalized Detector Counts")


def read_sed(path, encoding="latin-1"):
    """Parse one .sed file. Raises ValueError with a readable message on bad input.

    ``encoding`` defaults to latin-1 because DARWin writes Windows-encoded degree and
    micro signs into the header. A UTF-8 BOM, if present, is stripped before decoding:
    decoded as latin-1 it would otherwise become three visible characters glued to the
    first header key, silently turning "Version" into something no lookup can find.
    """
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw.startswith(b"\xef\xbb\xbf"):
        raw = raw[3:]
    lines = raw.decode(encoding, errors="replace").splitlines()

    header, marker = {}, None
    for i, line in enumerate(lines):
        if _DATA.match(line):
            marker = i
            break
        if ":" in line:
            k, _, v = line.partition(":")
            header[k.strip()] = v.strip()

    if marker is None:
        first = lines[0][:70] if lines else "<empty file>"
        raise ValueError(
            "%s has no 'Data:' line, so it is not a DARWin .sed export.\n"
            "First line was: %r" % (os.path.basename(path), first))

    body = [ln for ln in lines[marker + 1:] if ln.strip()]
    if len(body) < 2:
        raise ValueError("%s: no data rows after 'Data:'" % os.path.basename(path))

    labels = [c.strip() for c in body[0].split("\t")]
    if len(labels) < 2:
        labels = [c.strip() for c in re.split(r"\s{2,}|\t", body[0].strip())]

    rows = []
    for ln in body[1:]:
        parts = ln.split("\t") if "\t" in ln else ln.split()
        if len(parts) != len(labels):
            continue
        try:
            rows.append([float(p) for p in parts])
        except ValueError:
            continue
    if not rows:
        raise ValueError(
            "%s: header %s parsed but no numeric rows followed"
            % (os.path.basename(path), labels))

    raw = {}
    for j, lab in enumerate(labels):
        raw[lab] = [r[j] for r in rows]

    columns, scale = {}, None
    for canon, aliases in COLUMN_ALIASES.items():
        for lab in labels:
            if lab in aliases:
                columns[canon] = raw[lab]
                if canon == "reflectance":
                    scale = 100.0 if "%" in lab else 1.0
                break

    if "wavelength" not in columns:
        raise ValueError(
            "%s: no 'Wvl' wavelength column. Columns found: %s"
            % (os.path.basename(path), labels))

    order = sorted(range(len(columns["wavelength"])),
                   key=lambda i: columns["wavelength"][i])
    if order != list(range(len(order))):
        columns = {k: [v[i] for i in order] for k, v in columns.items()}
        raw = {k: [v[i] for i in order] for k, v in raw.items()}

    return SedSpectrum(path, header, columns, raw, scale)


def read_folder(folder, pattern=".sed"):
    """Read every .sed in a folder, sorted by name. Returns (spectra, errors)."""
    out, errs = [], []
    for fn in sorted(os.listdir(folder)):
        if not fn.lower().endswith(pattern):
            continue
        p = os.path.join(folder, fn)
        try:
            out.append(read_sed(p))
        except Exception as exc:
            errs.append((fn, "%s: %s" % (type(exc).__name__, exc)))
    return out, errs


def guess_role(spec):
    """Guess water / sky / panel from the filename and Comment field.

    Only a convenience for pre-selecting the GUI dropdowns. The operator always
    confirms, because a mislabelled sky scan silently produces a wrong R_rs.
    """
    text = ("%s %s" % (spec.name, spec.comment)).lower()
    for role, keys in (("panel", ("panel", "plaque", "spectralon", "white", "ref")),
                       ("sky", ("sky", "skye", "zenith")),
                       ("water", ("water", "sea", "target", "tgt", "station", "lake"))):
        for k in keys:
            if k in text:
                return role
    return "unassigned"


def _coord(value):
    if not value or value.strip().lower() in ("n/a", "", "none"):
        return None
    try:
        return float(re.sub(r"[^0-9.eE+-]", "", value.strip()))
    except ValueError:
        return None
