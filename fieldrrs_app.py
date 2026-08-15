"""Entry point for the frozen executable.

PyInstaller needs a plain script to freeze, not a ``python -m package`` invocation, so
this exists purely as the target of the build. It is the same GUI you get from
``python -m fieldrrs``.
"""

import multiprocessing
import sys


def main():
    # Harmless on a single-process app, and required if the frozen exe is ever run on a
    # machine where a child process gets spawned: without it a frozen Windows build can
    # re-launch itself in a loop.
    multiprocessing.freeze_support()

    try:
        from fieldrrs.gui import main as gui_main
    except Exception as exc:                       # pragma: no cover - startup guard
        # A frozen windowed build has no console, so an import failure would otherwise
        # be a silent no-op: the user double-clicks and nothing happens at all.
        try:
            import tkinter as tk
            from tkinter import messagebox
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror(
                "fieldrrs failed to start",
                "%s: %s\n\nThis usually means the build is incomplete. Rebuild with "
                "build_exe.bat, or run the source version with run_gui.bat."
                % (type(exc).__name__, exc))
        except Exception:
            print("fieldrrs failed to start: %s: %s" % (type(exc).__name__, exc))
        return 1

    gui_main()
    return 0


if __name__ == "__main__":
    sys.exit(main())
