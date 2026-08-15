"""fieldrrs: field spectroradiometer scans to remote-sensing reflectance.

Zero dependencies. Standard library only, so it runs on any Python 3.8+ install,
including a bare python.org install on a Windows field tablet with no internet.

    from fieldrrs import read_sed, rrs_from_sed
    water = read_sed("water_001.sed")
    sky   = read_sed("sky_001.sed")
    res   = rrs_from_sed(water, sky, rho=0.028, panel_reflectance=0.99)
    print(res.notes)          # read these before using the numbers

GUI:  python -m fieldrrs      (or double-click run_gui.bat on Windows)

The measurement protocol these equations assume is in FIELD_PROTOCOL.md. Read it before
collecting data; no amount of processing rescues the wrong three scans.
"""

from .resample import (
    SATELLITE_BANDS,
    bin_spectrum,
    gaussian_resample,
    write_batch_csv,
    write_rrs_csv,
)
from .rrs import (
    DEFAULT_PANEL_REFLECTANCE,
    RHO_MOBLEY1999,
    RrsResult,
    average_results,
    residual_correction,
    overcast_notes,
    rho_advice,
    rrs_from_sed,
    rrs_three_scan,
)
from .sed import SedSpectrum, guess_role, read_folder, read_sed

__version__ = "1.0.0"

__all__ = [
    "read_sed", "read_folder", "SedSpectrum", "guess_role",
    "rrs_from_sed", "rrs_three_scan", "residual_correction", "rho_advice",
    "average_results", "RrsResult", "RHO_MOBLEY1999", "DEFAULT_PANEL_REFLECTANCE",
    "overcast_notes",
    "bin_spectrum", "gaussian_resample", "write_rrs_csv", "write_batch_csv",
    "SATELLITE_BANDS", "__version__",
]
