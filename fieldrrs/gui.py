"""Field GUI: .sed scans in, R_rs out. Tkinter only, no third-party packages.

Designed for a Windows tablet in daylight: large fonts, large hit targets, warnings
that cannot be dismissed, and no step that silently guesses on the operator's behalf.

Launch:  python -m fieldrrs      or double-click run_gui.bat
"""

from __future__ import annotations

import os
import traceback

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from . import resample as rs
from .rrs import (
    DEFAULT_PANEL_REFLECTANCE,
    RHO_MOBLEY1999,
    average_results,
    rho_advice,
    rrs_from_sed,
)
from .sed import guess_role, read_folder, read_sed
from .solar import (
    ALT_RELATIVE_AZIMUTH,
    DEFAULT_RELATIVE_AZIMUTH,
    DEFAULT_VIEW_ZENITH,
    declination_from_sun_sighting,
    local_to_utc_hours,
    pointing,
)

BIG = ("Segoe UI", 12)
BIGB = ("Segoe UI", 12, "bold")
HUGE = ("Segoe UI", 14, "bold")
MONO = ("Consolas", 10)

ROLES = ("water", "sky", "panel", "skip")
ROLE_COLOR = {"water": "#0b6", "sky": "#08c", "panel": "#c60", "skip": "#999",
              "unassigned": "#333"}


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("fieldrrs - above-water Rrs from Spectral Evolution scans")
        self.geometry("1400x900")
        self.scans = []          # list of SedSpectrum  (batch mode)
        self.roles = []          # parallel list of role strings
        self.results = []        # list of (station_name, RrsResult)
        self.folder = None
        self.slots = {"water": None, "sky": None, "panel": None}   # single-station mode
        self.last_pointing = None
        self.last_footprint = None
        self._build()

    # ---------------------------------------------------------------- layout
    def _build(self):
        root = ttk.Frame(self, padding=6)
        root.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(root)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 8))
        right = ttk.Frame(root)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # --- step 1a: the three scans, loaded explicitly one at a time
        ttk.Label(left, text="1.  LOAD THE THREE SCANS", font=HUGE).pack(anchor="w")
        slots = ttk.Frame(left)
        slots.pack(fill=tk.X, pady=(2, 4))
        self.slot_labels = {}
        for i, (role, hint) in enumerate((
                ("water", "40 deg from NADIR, 135 deg from sun"),
                ("sky", "40 deg from ZENITH, same bearing as water"),
                ("panel", "level white panel  (optional)"))):
            tk.Button(slots, text="LOAD %s" % role.upper(), font=BIGB, width=12, height=2,
                      bg=ROLE_COLOR[role], fg="white",
                      command=lambda r=role: self.load_slot(r)).grid(
                          row=i, column=0, sticky="w", pady=2, padx=(0, 6))
            lab = tk.Label(slots, text="(not loaded)", font=("Consolas", 10),
                           fg="#a00", anchor="w", justify="left")
            lab.grid(row=i, column=1, sticky="w")
            self.slot_labels[role] = lab
            tk.Label(slots, text=hint, font=("Segoe UI", 8), fg="#555",
                     anchor="w").grid(row=i, column=2, sticky="w", padx=(8, 0))
        tk.Button(slots, text="clear all three", font=("Segoe UI", 9),
                  command=self.clear_slots).grid(row=3, column=0, columnspan=2,
                                                 sticky="w", pady=(2, 0))
        ttk.Label(left, text="PANEL is optional: leave it empty and the panel radiance "
                             "is read from the water file's own Rad. (Ref.) column\n"
                             "(the DARWin reference-scan workflow).",
                  font=("Segoe UI", 8), foreground="#444").pack(anchor="w")

        ttk.Separator(left).pack(fill=tk.X, pady=6)

        # --- step 1b: or a whole folder at once, for many stations
        b = tk.Button(left, text="OR:  OPEN WHOLE FOLDER (batch)", font=BIGB, height=2,
                      bg="#08c", fg="white", command=self.open_folder)
        b.pack(fill=tk.X, pady=(0, 4))

        ttk.Label(left, text="2.  Assign each scan  (batch mode only)",
                  font=BIGB).pack(anchor="w")
        ttk.Label(left, text="Select rows, then tap a role button.",
                  font=("Segoe UI", 9)).pack(anchor="w")

        cols = ("role", "file", "comment", "time")
        self.tree = ttk.Treeview(left, columns=cols, show="headings", height=16,
                                 selectmode="extended")
        for c, w in zip(cols, (80, 210, 150, 130)):
            self.tree.heading(c, text=c.upper())
            self.tree.column(c, width=w, anchor="w")
        self.tree.pack(fill=tk.X)
        style = ttk.Style(self)
        style.configure("Treeview", font=("Segoe UI", 10), rowheight=26)

        rowb = ttk.Frame(left)
        rowb.pack(fill=tk.X, pady=4)
        for role in ROLES:
            tk.Button(rowb, text=role.upper(), font=BIGB, width=7, height=2,
                      bg=ROLE_COLOR[role], fg="white",
                      command=lambda r=role: self.assign(r)).pack(side=tk.LEFT, padx=2)

        # --- step 3: settings
        ttk.Label(left, text="3.  Settings", font=BIGB).pack(anchor="w", pady=(8, 0))
        g = ttk.Frame(left)
        g.pack(fill=tk.X)
        self.panel_r = self._entry(g, "Panel reflectance", str(DEFAULT_PANEL_REFLECTANCE), 0)
        self.wind = self._entry(g, "Wind speed (m/s)", "", 1)
        self.rho = self._entry(g, "rho (sky glint)", str(RHO_MOBLEY1999), 2)
        tk.Button(g, text="check rho vs wind", font=("Segoe UI", 9),
                  command=self.check_rho).grid(row=3, column=0, columnspan=2,
                                               sticky="ew", pady=2)
        self.residual = self._combo(
            g, "Residual glint", ["none", "nir_zero", "nir_similarity"], 4)
        self.source = self._combo(g, "Use columns", ["radiance", "reflectance"], 5)

        # --- geometry + where is the sun
        ttk.Label(left, text="Geometry  (defaults = Mobley 1999)",
                  font=BIGB).pack(anchor="w", pady=(8, 0))
        gg = ttk.Frame(left)
        gg.pack(fill=tk.X)
        self.view_zen = self._entry(gg, "View zenith (from nadir)",
                                    str(DEFAULT_VIEW_ZENITH), 0)
        self.rel_az = self._combo(gg, "Rel. azimuth from sun",
                                  [str(DEFAULT_RELATIVE_AZIMUTH),
                                   str(ALT_RELATIVE_AZIMUTH)], 1)
        self.lat = self._entry(gg, "Latitude (N +)", "", 2)
        self.lon = self._entry(gg, "Longitude (E +)", "", 3)
        self.date = self._entry(gg, "Date  YYYY-MM-DD", "", 4)
        self.clock = self._entry(gg, "Local time  HH:MM", "", 5)
        self.utc_off = self._entry(gg, "UTC offset (h)", "0", 6)
        self.mag_dec = self._entry(gg, "Mag. declination (deg)", "", 9)
        tk.Button(gg, text="WHERE IS THE SUN?  ->  bearings", font=BIGB, height=2,
                  bg="#c60", fg="white", command=self.where_is_the_sun).grid(
                      row=7, column=0, columnspan=2, sticky="ew", pady=4)
        tk.Button(gg, text="fill from the loaded water scan", font=("Segoe UI", 9),
                  command=self.fill_from_scan).grid(row=8, column=0, columnspan=2,
                                                    sticky="ew")

        # --- step 4: run
        tk.Button(left, text="4.  COMPUTE  Rrs", font=HUGE, height=2,
                  bg="#0a6", fg="white", command=self.compute).pack(fill=tk.X, pady=6)

        outb = ttk.Frame(left)
        outb.pack(fill=tk.X)
        tk.Button(outb, text="Save this spectrum", font=BIG, height=2,
                  command=self.save_one).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(outb, text="Save ALL (batch)", font=BIG, height=2,
                  command=self.save_all).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

        # --- right: plot + messages
        self.canvas = tk.Canvas(right, bg="white", height=470,
                                highlightthickness=1, highlightbackground="#999")
        self.canvas.pack(fill=tk.BOTH, expand=True)
        self.canvas.bind("<Configure>", lambda e: self.redraw())

        ttk.Label(right, text="Messages  (read these before trusting a number)",
                  font=BIGB).pack(anchor="w", pady=(6, 0))
        self.log = tk.Text(right, height=11, font=MONO, wrap=tk.WORD, bg="#fffdf3")
        self.log.pack(fill=tk.BOTH, expand=False)

        ttk.Label(right, text="Footprint  (what patch of water this scan averaged over)",
                  font=BIGB).pack(anchor="w", pady=(6, 0))
        self.fp_canvas = tk.Canvas(right, bg="white", height=190,
                                   highlightthickness=1, highlightbackground="#999")
        self.fp_canvas.pack(fill=tk.X)
        self.fp_canvas.bind("<Configure>", lambda e: self.draw_footprint())
        self.say("Ready. Open the folder of .sed scans from today's station.")
        self.say("Reminder: water = 40 deg from nadir, 135 deg from the sun; "
                 "sky = 40 deg from ZENITH at the same bearing; panel level and unshaded.")

        self.current = None

    def _entry(self, parent, label, default, row):
        ttk.Label(parent, text=label, font=BIG).grid(row=row, column=0, sticky="w", pady=2)
        v = tk.StringVar(value=default)
        tk.Entry(parent, textvariable=v, font=BIG, width=10).grid(row=row, column=1,
                                                                  sticky="e")
        return v

    def _combo(self, parent, label, values, row):
        ttk.Label(parent, text=label, font=BIG).grid(row=row, column=0, sticky="w", pady=2)
        v = tk.StringVar(value=values[0])
        ttk.Combobox(parent, textvariable=v, values=values, font=BIG, width=13,
                     state="readonly").grid(row=row, column=1, sticky="e")
        return v

    # --------------------------------------------------------------- actions
    def load_slot(self, role):
        """Load one named scan. The operator says which is which; nothing is guessed."""
        p = filedialog.askopenfilename(
            title="Select the %s scan" % role.upper(),
            initialdir=self.folder or os.getcwd(),
            filetypes=[("Spectral Evolution scan", "*.sed"), ("All files", "*.*")])
        if not p:
            return
        try:
            spec = read_sed(p)
        except Exception as exc:
            return self._err(exc)

        guessed = guess_role(spec)
        if guessed != "unassigned" and guessed != role:
            if not messagebox.askyesno(
                    "Name does not match",
                    "You are loading this as the %s scan, but its name/comment says "
                    "'%s':\n\n%s\ncomment: %s\n\nLoad it as %s anyway?"
                    % (role.upper(), guessed, spec.name, spec.comment, role.upper())):
                return

        self.slots[role] = spec
        self.folder = os.path.dirname(p)
        self.slot_labels[role].config(
            text="%s   (%d bands %.0f-%.0f nm)"
                 % (spec.name, len(spec.wavelength), spec.wavelength[0],
                    spec.wavelength[-1]),
            fg="#060")
        self.say("Loaded %s scan: %s   %s" % (role.upper(), spec.name, spec.when))
        self._report_instrument_metadata(spec, role)
        if role == "water" and self.slots["panel"] is None:
            if spec.has("rad_ref"):
                self.say("   panel radiance will come from this file's Rad. (Ref.) "
                         "column unless you load a separate PANEL scan.")
            else:
                self.say("   ! this file has no Rad. (Ref.) column, so a separate "
                         "PANEL scan IS required.")

    def fill_from_scan(self):
        """Take lat/lon/date/time from the loaded water scan's own header if it has them."""
        spec = self.slots["water"] or (self._picked("water") or [None])[0]
        if spec is None:
            return messagebox.showinfo("No water scan", "Load a WATER scan first.")
        got = []
        if spec.latitude is not None:
            self.lat.set("%.5f" % spec.latitude); got.append("lat")
        if spec.longitude is not None:
            self.lon.set("%.5f" % spec.longitude); got.append("lon")
        d = spec.header.get("Date", "").split(",")[0].strip()
        t = spec.header.get("Time", "").split(",")[0].strip()
        if d:
            parts = d.replace("-", "/").split("/")
            if len(parts) == 3:
                # DARWin writes MM/DD/YYYY
                mm, dd, yy = parts
                if len(parts[0]) == 4:
                    yy, mm, dd = parts
                try:
                    self.date.set("%04d-%02d-%02d" % (int(yy), int(mm), int(dd)))
                    got.append("date")
                except ValueError:
                    pass
        if t and ":" in t:
            self.clock.set(":".join(t.split(":")[:2])); got.append("time")
        self.say("Filled %s from %s. CHECK the UTC offset: the instrument clock is "
                 "local time and the header does not say which zone."
                 % (", ".join(got) or "nothing", spec.name))

    def where_is_the_sun(self):
        try:
            lat = float(self.lat.get()); lon = float(self.lon.get())
            y, mo, d = [int(x) for x in self.date.get().strip().split("-")]
            hh, mm = [int(x) for x in self.clock.get().strip().split(":")]
            off = float(self.utc_off.get() or 0.0)
            vz = float(self.view_zen.get()); ra = float(self.rel_az.get())
        except Exception:
            return messagebox.showerror(
                "Geometry inputs",
                "Need latitude, longitude, date as YYYY-MM-DD, local time as HH:MM, "
                "and the UTC offset in hours (e.g. -4 for US Eastern daylight time).")

        hour_utc, shift = local_to_utc_hours(hh, mm, off)
        import datetime as _dt
        dd = _dt.date(y, mo, d) + _dt.timedelta(days=shift)
        p = pointing(lat, lon, dd.year, dd.month, dd.day, hour_utc,
                     view_zenith=vz, relative_azimuth=ra)
        self.last_pointing = p

        self.say("")
        self.say("=== POINTING  (%04d-%02d-%02d %02d:%02d local, %.2f h UTC) ==="
                 % (dd.year, dd.month, dd.day, hh, mm, hour_utc))
        dec = None
        try:
            dec = float(self.mag_dec.get()) if self.mag_dec.get().strip() else None
        except ValueError:
            self.say("   (magnetic declination not a number, ignoring)")
        for line in p.describe(declination=dec):
            self.say("   " + line)
        if dec is None:
            self.say("   TIP: point your phone at the sun, read the magnetic bearing,")
            self.say("        then set Mag. declination = %.0f minus that reading."
                     % p.sun.azimuth)
        self.say("   " + p.sun.advice())
        if not p.sun.usable:
            messagebox.showwarning("Solar zenith", p.sun.advice())

    def _report_instrument_metadata(self, spec, role):
        """Surface what the instrument itself logged: attitude, solar angle, footprint."""
        uf = spec.user_fields
        if uf:
            self.say("   instrument logged: " +
                     ", ".join("%s %s" % (k, v) for k, v in uf.items()))
        if spec.solar_elevation_deg is not None:
            self.say("   Solar Angle %.2f deg is solar ELEVATION (verified against an "
                     "independent calculation)." % spec.solar_elevation_deg)
        # The instrument clock is set by hand and can be wrong; GPS time is UTC.
        gps = spec.gps_time
        if gps is not None:
            self.say("   GPS time %.2f h UTC -- USE THIS for solar geometry, not the "
                     "instrument clock." % gps)
        try:
            vz = float(self.view_zen.get())
        except Exception:
            vz = 40.0
        fp = spec.footprint(vz)
        if fp is not None and role == "water":
            self.last_footprint = fp
            self.say("   FOOTPRINT at %.0f deg: %.2f x %.2f m ellipse, %.2f m2, "
                     "sensor %.2f m above the surface, spot %.2f m out."
                     % (vz, fp["spot_across_m"], fp["spot_along_m"], fp["area_m2"],
                        fp["height_above_surface_m"], fp["horizontal_offset_m"]))
            self.say("   (an OLCI 300 m pixel is %.0fx that area; MODIS 1 km, %.0fx)"
                     % (300*300/fp["area_m2"], 1000*1000/fp["area_m2"]))
            self.draw_footprint()

    def draw_footprint(self):
        """To-scale plan view of the water patch this scan actually averaged over."""
        fp = getattr(self, "last_footprint", None)
        c = self.fp_canvas
        c.delete("all")
        W = c.winfo_width() or 360
        H = c.winfo_height() or 190
        if not fp:
            c.create_text(W/2, H/2, text="footprint appears when a WATER scan is loaded",
                          font=("Segoe UI", 9), fill="#888")
            return
        a, b = fp["spot_along_m"], fp["spot_across_m"]      # along-view, across-view
        span = max(a, b, 1.0) * 1.45
        sc = min(W - 60, H - 46) / span
        cx, cy = W/2, H/2 + 6

        # 1 m scale bar
        c.create_line(14, H-14, 14 + sc, H-14, fill="#333", width=3)
        c.create_text(14 + sc/2, H-26, text="1 m", font=("Consolas", 9))

        c.create_oval(cx - a*sc/2, cy - b*sc/2, cx + a*sc/2, cy + b*sc/2,
                      outline="#0b6", width=3, fill="#d9f2e4")
        c.create_text(cx, cy - 4, text="%.2f x %.2f m" % (a, b),
                      font=("Segoe UI", 10, "bold"), fill="#064")
        c.create_text(cx, cy + 14, text="%.2f m\u00b2" % fp["area_m2"],
                      font=("Segoe UI", 9), fill="#064")
        c.create_text(W/2, 14,
                      text="water patch sampled  (view %.0f deg, range %.2f m, FOV %.0f deg)"
                           % (fp["view_zenith_deg"], fp["range_m"], fp["fov_deg"]),
                      font=("Segoe UI", 9))
        # instrument position, projected
        ix = cx - fp["horizontal_offset_m"] * sc
        if ix > 8:
            c.create_line(ix, cy, cx, cy, fill="#888", dash=(3, 3))
            c.create_oval(ix-4, cy-4, ix+4, cy+4, fill="#333")
            c.create_text(ix, cy - 16, text="sensor\n(%.1f m up)" % fp["height_above_surface_m"],
                          font=("Consolas", 8), fill="#333")

    def clear_slots(self):
        for role in self.slots:
            self.slots[role] = None
            self.slot_labels[role].config(text="(not loaded)", fg="#a00")
        self.say("Cleared the three scan slots.")

    def open_folder(self):
        d = filedialog.askdirectory(title="Folder containing today's .sed scans")
        if not d:
            return
        try:
            scans, errs = read_folder(d)
        except Exception as exc:
            return self._err(exc)
        if not scans:
            messagebox.showwarning(
                "No scans", "No readable .sed files in:\n%s\n\n%s" %
                (d, "\n".join("%s: %s" % e for e in errs) if errs else ""))
            return
        self.folder, self.scans = d, scans
        self.roles = [guess_role(s) for s in scans]
        self.refresh_tree()
        self.say("Loaded %d scans from %s" % (len(scans), d))
        for fn, msg in errs:
            self.say("SKIPPED %s -> %s" % (fn, msg))
        counts = {r: self.roles.count(r) for r in set(self.roles)}
        self.say("Auto-guessed roles from filename/comment: %s" % counts)
        self.say("CHECK THESE. A sky scan mislabelled as water gives a wrong answer "
                 "with no error message.")

    def refresh_tree(self):
        self.tree.delete(*self.tree.get_children())
        for i, (s, role) in enumerate(zip(self.scans, self.roles)):
            self.tree.insert("", "end", iid=str(i),
                             values=(role, s.name, s.comment[:28], s.when))

    def assign(self, role):
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo("Select first", "Select one or more rows, then tap a role.")
            return
        for iid in sel:
            self.roles[int(iid)] = role
        self.refresh_tree()
        for iid in sel:
            self.tree.selection_add(iid)

    def check_rho(self):
        txt = self.wind.get().strip()
        try:
            wind = float(txt) if txt else None
        except ValueError:
            return messagebox.showerror("Wind", "Wind speed must be a number (m/s).")
        val, msg = rho_advice(wind)
        self.say(msg)
        if val is None:
            messagebox.showwarning("rho outside its validity range", msg)
        else:
            self.rho.set(str(val))

    def _picked(self, role):
        return [s for s, r in zip(self.scans, self.roles) if r == role]

    def compute(self):
        try:
            # The three explicit slots win when they are filled; otherwise fall back to
            # whatever the batch table has been assigned.
            if self.slots["water"] is not None or self.slots["sky"] is not None:
                waters = [self.slots["water"]] if self.slots["water"] else []
                skies = [self.slots["sky"]] if self.slots["sky"] else []
                panels = [self.slots["panel"]] if self.slots["panel"] else []
                mode = "three loaded scans"
            else:
                waters, skies, panels = (self._picked("water"), self._picked("sky"),
                                         self._picked("panel"))
                mode = "folder / batch"

            if not waters:
                return messagebox.showwarning(
                    "No water scan",
                    "Load a WATER scan (button 1), or assign one in the batch table.")
            if not skies:
                return messagebox.showwarning(
                    "No sky scan",
                    "Load a SKY scan.\n\nThe sky scan is NOT optional. Without it the "
                    "reflected skylight stays in the signal and R_rs comes out far too "
                    "high in the blue.\n\nIt is the sky at 40 deg from ZENITH on the "
                    "SAME compass bearing as your water view.")
            sky = skies[0]
            panel = panels[0] if panels else None
            self.say("Computing from: %s" % mode)
            if panel is None and self.source.get() == "radiance":
                self.say("No panel scan assigned; using the water file's own "
                         "'Rad. (Ref.)' column as the panel radiance (the DARWin "
                         "reference-scan workflow).")

            pr = float(self.panel_r.get())
            rho = float(self.rho.get())
            resid = self.residual.get()
            src = self.source.get()

            self.results = []
            self.log.delete("1.0", tk.END)
            geom = self._geometry_meta()
            for w in waters:
                res = rrs_from_sed(w, sky, panel, panel_reflectance=pr, rho=rho,
                                   residual=resid, source=src)
                res.meta.update(geom)
                res.meta["panel_reflectance"] = pr
                self.results.append((os.path.splitext(w.name)[0], res))
                self.say("--- %s" % w.name)
                for n in res.notes:
                    self.say("   ! " + n)
                for band in (443, 490, 555, 670):
                    wlv, v = res.value_at(band)
                    self.say("   Rrs(%.0f) = %.5f sr^-1" % (wlv, v))

            if len(self.results) > 1:
                try:
                    avg = average_results([r for _, r in self.results])
                    self.say("Replicate mean of %d water scans computed." % len(self.results))
                    self.results.append(("MEAN_of_%d" % (len(self.results)), avg))
                except Exception as exc:
                    self.say("Could not average replicates: %s" % exc)

            self.current = self.results[-1][1]
            self.redraw()
            self.say("Done. %d spectra. Check the warnings above." % len(self.results))
        except Exception as exc:
            self._err(exc)

    def _geometry_meta(self):
        """Everything about how the measurement was aimed, written into every CSV.

        An R_rs without its viewing geometry and wind speed cannot be reprocessed later,
        because rho is conditional on both and there is no way to recover them.
        """
        meta = {
            "view_zenith_from_nadir_deg": self.view_zen.get(),
            "relative_azimuth_from_sun_deg": self.rel_az.get(),
            "wind_speed_ms": self.wind.get() or "NOT RECORDED",
        }
        fp = getattr(self, "last_footprint", None)
        if fp is not None:
            meta.update({
                "footprint_across_m": "%.3f" % fp["spot_across_m"],
                "footprint_along_m": "%.3f" % fp["spot_along_m"],
                "footprint_area_m2": "%.3f" % fp["area_m2"],
                "sensor_height_m": "%.3f" % fp["height_above_surface_m"],
                "range_m": "%.3f" % fp["range_m"],
                "fov_deg": "%.1f" % fp["fov_deg"],
            })
        p = self.last_pointing
        if p is not None:
            meta.update({
                "solar_zenith_deg": "%.2f" % p.sun.zenith,
                "solar_azimuth_deg": "%.2f" % p.sun.azimuth,
                "solar_elevation_deg": "%.2f" % p.sun.elevation,
                "view_bearing_options_deg": "%.0f or %.0f"
                                            % (p.bearing_ccw, p.bearing_cw),
                "solar_zenith_in_usable_window": str(p.sun.usable),
            })
        else:
            meta["solar_geometry"] = ("NOT COMPUTED - press WHERE IS THE SUN before "
                                      "computing to record it")
        return meta

    def save_one(self):
        if not self.current:
            return messagebox.showwarning("Nothing yet", "Compute first.")
        p = filedialog.asksaveasfilename(
            defaultextension=".csv", initialdir=self.folder,
            initialfile="rrs_spectrum.csv", filetypes=[("CSV", "*.csv")])
        if p:
            rs.write_rrs_csv(p, self.current)
            self.say("Wrote %s" % p)

    def save_all(self):
        if not self.results:
            return messagebox.showwarning("Nothing yet", "Compute first.")
        d = filedialog.askdirectory(title="Where to write the results",
                                    initialdir=self.folder)
        if not d:
            return
        for name, res in self.results:
            rs.write_rrs_csv(os.path.join(d, "rrs_%s.csv" % name), res)
        batch = os.path.join(d, "rrs_all_stations.csv")
        rs.write_batch_csv(batch, [(n, r) for n, r in self.results])
        self.say("Wrote %d spectra + %s" % (len(self.results), batch))
        messagebox.showinfo("Saved", "Wrote %d files to\n%s" % (len(self.results) + 1, d))

    # ------------------------------------------------------------- plotting
    def redraw(self):
        c = self.canvas
        c.delete("all")
        W = c.winfo_width() or 800
        H = c.winfo_height() or 470
        L, R, T, B = 80, 20, 24, 46
        if not self.current:
            c.create_text(W / 2, H / 2, text="R_rs will be plotted here",
                          font=BIG, fill="#888")
            return

        series = [(n, r) for n, r in self.results]
        xs, ys = [], []
        for _, r in series:
            for w, v in zip(r.wavelength, r.rrs):
                if v == v and 350 <= w <= 1000:
                    xs.append(w); ys.append(v)
        if not xs:
            c.create_text(W / 2, H / 2, text="no finite values to plot", font=BIG)
            return
        x0, x1 = min(xs), max(xs)
        y0, y1 = min(ys), max(ys)
        if y1 <= y0:
            y1 = y0 + 1e-6
        pad = 0.08 * (y1 - y0)
        y0, y1 = y0 - pad, y1 + pad

        def px(w):
            return L + (w - x0) / (x1 - x0) * (W - L - R)

        def py(v):
            return H - B - (v - y0) / (y1 - y0) * (H - T - B)

        c.create_rectangle(L, T, W - R, H - B, outline="#666")
        for k in range(6):
            v = y0 + k * (y1 - y0) / 5.0
            yy = py(v)
            c.create_line(L, yy, W - R, yy, fill="#eee")
            c.create_text(L - 6, yy, text="%.4f" % v, anchor="e", font=("Consolas", 9))
        for w in range(int(x0 // 100 * 100), int(x1) + 100, 100):
            if x0 <= w <= x1:
                c.create_line(px(w), T, px(w), H - B, fill="#eee")
                c.create_text(px(w), H - B + 14, text=str(w), font=("Consolas", 9))
        if y0 < 0 < y1:
            c.create_line(L, py(0), W - R, py(0), fill="#c00", dash=(4, 3))
            c.create_text(W - R - 4, py(0) - 8, text="Rrs = 0", anchor="e",
                          fill="#c00", font=("Consolas", 9))

        palette = ["#0a6", "#08c", "#c60", "#a0a", "#666", "#c00"]
        for i, (name, r) in enumerate(series):
            col = palette[i % len(palette)]
            pts = []
            for w, v in zip(r.wavelength, r.rrs):
                if v == v and x0 <= w <= x1:
                    pts.extend((px(w), py(v)))
            if len(pts) >= 4:
                c.create_line(*pts, fill=col, width=3 if "MEAN" in name else 2)
            c.create_text(W - R - 8, T + 14 + 16 * i, text=name[:26], anchor="e",
                          fill=col, font=("Consolas", 10))

        c.create_text((L + W - R) / 2, H - 8, text="wavelength (nm)", font=BIG)
        c.create_text(16, (T + H - B) / 2, text="Rrs (1/sr)", font=BIG, angle=90)

    # ----------------------------------------------------------------- misc
    def say(self, msg):
        self.log.insert(tk.END, msg + "\n")
        self.log.see(tk.END)

    def _err(self, exc):
        messagebox.showerror(type(exc).__name__, "%s\n\n%s" % (exc, traceback.format_exc()))
        self.say("ERROR %s: %s" % (type(exc).__name__, exc))


def main():
    App().mainloop()


if __name__ == "__main__":
    main()
