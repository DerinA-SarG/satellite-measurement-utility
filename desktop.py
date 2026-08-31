"""Desktop wrapper for Satellite Measurement Utility.

Serves the web assets from a loopback-only HTTP server and shows them in a
native window (WebView2 on Windows). Saving and opening files go through real
Windows file dialogs rather than browser downloads.

Run directly for development:      python desktop.py
Build a standalone exe:            python build_exe.py
"""

import datetime
import functools
import glob
import http.server
import os
import socketserver
import sys
import tempfile
import threading
import traceback

APP_NAME = "Satellite Measurement Utility"
BACKUP_KEEP = 40        # newest N of each format; older ones are pruned

# The page hands its exported text over as it goes, so that closing the window
# only has to write a file. Pulling it out of the webview at close time instead
# would mean calling into the UI thread while it is being torn down.
_stash = {"kml": None, "geojson": None, "shapes": 0}


def asset_dir():
    """Where index.html lives — the PyInstaller bundle, or this folder."""
    return getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))


def app_dir():
    """The folder the exe itself sits in, not the unpacked bundle."""
    if getattr(sys, "frozen", False):
        return os.path.dirname(os.path.abspath(sys.executable))
    return os.path.dirname(os.path.abspath(__file__))


def backup_dir():
    """Beside the exe if that is writable, otherwise under Documents. Returns
    None if neither is, in which case there is nowhere to back up to."""
    for cand in (os.path.join(app_dir(), "backups"),
                 os.path.join(os.path.expanduser("~"), "Documents", APP_NAME, "backups")):
        try:
            os.makedirs(cand, exist_ok=True)
            probe = os.path.join(cand, ".writable")
            with open(probe, "w", encoding="utf-8") as fh:
                fh.write("")
            os.remove(probe)
            return cand
        except OSError:
            continue
    return None


def prune(folder):
    """Keep the newest BACKUP_KEEP of each format."""
    for ext in ("kml", "geojson"):
        found = sorted(glob.glob(os.path.join(folder, "site-measurements-*." + ext)),
                       key=os.path.getmtime, reverse=True)
        for old in found[BACKUP_KEEP:]:
            try:
                os.remove(old)
            except OSError:
                pass


def write_backup():
    """Write what the page last handed over. The window keeps nothing between
    runs by design, so this is the copy that survives it."""
    if not _stash["shapes"] or not (_stash["kml"] or _stash["geojson"]):
        log("backup: nothing drawn, nothing to write")
        return None
    folder = backup_dir()
    if not folder:
        log("backup: nowhere writable to put it")
        return None
    stamp = datetime.datetime.now().strftime("%Y-%m-%d-%H%M%S")
    written = []
    for ext in ("kml", "geojson"):
        text = _stash[ext]
        if not text:
            continue
        path = os.path.join(folder, "site-measurements-{}.{}".format(stamp, ext))
        try:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(text)
            written.append(path)
        except OSError as e:
            log("backup {} failed: {}".format(ext, e))
    for path in written:
        log("backup -> " + path)
    prune(folder)
    return folder


def log_path():
    return os.path.join(tempfile.gettempdir(), "satellite-measurement-utility.log")


def log(msg):
    try:
        with open(log_path(), "a", encoding="utf-8") as fh:
            fh.write(msg.rstrip() + "\n")
    except OSError:
        pass


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def end_headers(self):
        # never let the webview cache a stale build
        self.send_header("Cache-Control", "no-store")
        super().end_headers()


def start_server(directory):
    """Serve `directory` on a free loopback port. Returns the port."""
    handler = functools.partial(QuietHandler, directory=directory)
    socketserver.TCPServer.allow_reuse_address = True
    httpd = socketserver.TCPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    return httpd.server_address[1]


class Api:
    """Exposed to the page as window.pywebview.api.

    Keep every public attribute a plain method: pywebview walks this object to
    build the JS bridge, so holding a reference to the window here sends it
    recursing through the WebView2 COM objects.
    """

    def stash(self, kml, geojson, shapes):
        """The page's current work, handed over whenever it autosaves. Held in
        memory only; write_backup puts it on disk as the window closes."""
        _stash["kml"] = kml
        _stash["geojson"] = geojson
        _stash["shapes"] = int(shapes or 0)
        return True

    def backup_folder(self):
        """Where a backup would go, so the page can say so. None if nowhere."""
        return backup_dir()

    def _dialog(self, kind, **kw):
        import webview
        result = webview.windows[0].create_file_dialog(kind, **kw)
        if not result:
            return None
        return result[0] if isinstance(result, (list, tuple)) else result

    def save_file(self, filename, content):
        """Native Save As. Returns the path written, or None if cancelled."""
        import webview
        try:
            ext = os.path.splitext(filename)[1].lower()
            types = (("KML file (*.kml)",) if ext == ".kml" else ("GeoJSON file (*.geojson)",)) \
                + ("All files (*.*)",)
            path = self._dialog(webview.SAVE_DIALOG, save_filename=filename, file_types=types)
            if not path:
                return None
            if not os.path.splitext(path)[1]:
                path += ext
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(content)
            return path
        except Exception:
            log(traceback.format_exc())
            return None

    def save_image(self, filename, b64):
        """Native Save As for the capture. Base64 in, because the bridge to the
        page carries text; returns the path written, or None if cancelled."""
        import base64
        import webview
        try:
            path = self._dialog(
                webview.SAVE_DIALOG,
                save_filename=filename,
                file_types=("PNG image (*.png)", "All files (*.*)"),
            )
            if not path:
                return None
            if not os.path.splitext(path)[1]:
                path += ".png"
            with open(path, "wb") as fh:
                fh.write(base64.b64decode(b64))
            return path
        except Exception:
            log(traceback.format_exc())
            return None

    def open_file(self):
        """Native Open. Returns {name, text} or None if cancelled."""
        import webview
        try:
            path = self._dialog(
                webview.OPEN_DIALOG,
                allow_multiple=False,
                file_types=("Map files (*.kml;*.geojson;*.json)", "All files (*.*)"),
            )
            if not path:
                return None
            with open(path, "r", encoding="utf-8") as fh:
                return {"name": os.path.basename(path), "text": fh.read()}
        except Exception:
            log(traceback.format_exc())
            return None


def main():
    try:
        import webview
    except ImportError:
        print("pywebview is not installed:  python -m pip install pywebview")
        return 1

    root = asset_dir()
    if not os.path.exists(os.path.join(root, "index.html")):
        log(f"index.html missing from {root}")
        print(f"index.html not found in {root}")
        return 1

    port = start_server(root)
    log(f"{APP_NAME} serving {root} on 127.0.0.1:{port}")

    # let a test harness find the port
    try:
        with open(os.path.join(tempfile.gettempdir(), "smu-port"), "w") as fh:
            fh.write(str(port))
    except OSError:
        pass

    window = webview.create_window(
        APP_NAME,
        f"http://127.0.0.1:{port}/",
        js_api=Api(),
        width=1360,
        height=880,
        min_size=(900, 620),
        text_select=False,
    )

    # Nothing survives the window otherwise -- it runs private, and the server
    # takes a new port each launch, so the page's own storage is a fresh one
    # every time. A copy in both formats on the way out is the safety net.
    def on_closing():
        log("closing: writing backup")
        try:
            folder = write_backup()
            if folder:
                log("backup folder: " + folder)
        except Exception:
            log(traceback.format_exc())
        return True     # never block the close over a backup

    window.events.closing += on_closing
    webview.start()
    return 0


if __name__ == "__main__":
    sys.exit(main())
