"""Build Satellite Measurement Utility into a single standalone .exe.

    python build_exe.py

Produces dist/SatelliteMeasurementUtility.exe — one file, no installer, no Python
needed on the machine that runs it. Needs pyinstaller, pywebview and Pillow
installed here to build.
"""

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NAME = "SatelliteMeasurementUtility"
ICON = os.path.join(HERE, "icon.ico")


def make_icon():
    """A roof outline with a measured corner, drawn at each icon size."""
    from PIL import Image, ImageDraw

    imgs = []
    for size in (16, 24, 32, 48, 64, 128, 256):
        img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
        d = ImageDraw.Draw(img)
        s = size / 32.0
        r = max(2, int(4 * s))
        d.rounded_rectangle([0, 0, size - 1, size - 1], radius=r, fill=(23, 30, 39, 255))

        pad = 6 * s
        box = [pad, pad + 2 * s, size - pad, size - pad - 2 * s]
        d.rectangle(box, fill=(56, 189, 248, 60),
                    outline=(56, 189, 248, 255), width=max(1, int(2 * s)))

        # corner handles
        h = max(1, int(2.2 * s))
        for x in (box[0], box[2]):
            for y in (box[1], box[3]):
                d.ellipse([x - h, y - h, x + h, y + h], fill=(255, 255, 255, 255))
        imgs.append(img)

    imgs[-1].save(ICON, format="ICO", sizes=[(i.width, i.height) for i in imgs])
    print(f"icon    -> {ICON}")


def main():
    os.chdir(HERE)
    make_icon()

    sep = ";" if os.name == "nt" else ":"
    data = ["index.html", "app.css", "app.js", "vendor"]
    missing = [d for d in data if not os.path.exists(os.path.join(HERE, d))]
    if missing:
        sys.exit(f"missing asset(s): {', '.join(missing)}")

    cmd = [sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean",
           "--onefile", "--windowed", "--name", NAME, "--icon", ICON]
    for d in data:
        cmd += ["--add-data", f"{d}{sep}{'vendor' if d == 'vendor' else '.'}"]
    cmd.append("desktop.py")

    print("build   -> " + " ".join(cmd[2:]))
    r = subprocess.run(cmd)
    if r.returncode != 0:
        sys.exit(r.returncode)

    for junk in ("build", f"{NAME}.spec"):
        p = os.path.join(HERE, junk)
        shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) else \
            (os.remove(p) if os.path.exists(p) else None)

    exe = os.path.join(HERE, "dist", f"{NAME}.exe")
    print(f"\ndone    -> {exe}  ({os.path.getsize(exe) / 1048576:.1f} MB)")


if __name__ == "__main__":
    main()
