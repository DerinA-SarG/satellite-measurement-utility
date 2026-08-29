"""Desktop wrapper for Satellite Measurement Utility.

Serves the web assets from a loopback-only HTTP server and shows them in a
native window (WebView2 on Windows). Saving and opening files go through real
Windows file dialogs rather than browser downloads.

Run directly for development:      python desktop.py
Build a standalone exe:            python build_exe.py
"""

import functools
import http.server
import os
import socketserver
import sys
import tempfile
import threading
import traceback

APP_NAME = "Satellite Measurement Utility"


def asset_dir():
    """Where index.html lives — the PyInstaller bundle, or this folder."""
    return getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))


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

    webview.create_window(
        APP_NAME,
        f"http://127.0.0.1:{port}/",
        js_api=Api(),
        width=1360,
        height=880,
        min_size=(900, 620),
        text_select=False,
    )
    webview.start()
    return 0


if __name__ == "__main__":
    sys.exit(main())
