# PyInstaller spec for fieldrrs.
#
#   pyinstaller --clean --noconfirm fieldrrs.spec
#
# Produces dist/fieldrrs.exe on Windows, a single self-contained file that needs no
# Python install. The package is pure standard library plus tkinter, so there is nothing
# to bundle beyond the interpreter itself and no hidden imports to chase.
#
# CONSOLE: set to True by default. A windowed build looks tidier, but if the exe fails
# to start in the field a windowed build shows you nothing at all. A console window that
# prints a traceback is worth more than a clean taskbar. Flip it if you prefer.

CONSOLE = True

a = Analysis(
    ['fieldrrs_app.py'],
    pathex=['.'],
    binaries=[],
    datas=[
        ('FIELD_PROTOCOL.md', '.'),
        ('README.md', '.'),
    ],
    hiddenimports=[
        'fieldrrs', 'fieldrrs.gui', 'fieldrrs.sed', 'fieldrrs.rrs',
        'fieldrrs.resample', 'fieldrrs.solar',
    ],
    hookspath=[],
    runtime_hooks=[],
    # Nothing scientific is bundled: the field package deliberately has no numpy/scipy/
    # matplotlib dependency, so excluding them keeps the exe small and the build fast.
    excludes=['numpy', 'scipy', 'matplotlib', 'pandas', 'PIL', 'pytest'],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='fieldrrs',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    runtime_tmpdir=None,
    console=CONSOLE,
    disable_windowed_traceback=False,
)
