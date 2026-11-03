"""Screenshot the site at phone width so the layout can be judged, not just measured."""
import os
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

ROOT = r"D:\site\docs"
OUT = r"D:\site\_mobile_shots"
PORT = 8898
SHOTS = [
    ("index.html", "home"),
    ("articles/systems.html", "systems"),
    ("articles/signal.html", "signal"),
]


def serve():
    handler = partial(SimpleHTTPRequestHandler, directory=ROOT)
    httpd = ThreadingHTTPServer(("127.0.0.1", PORT), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    return httpd


def main():
    from playwright.sync_api import sync_playwright

    os.makedirs(OUT, exist_ok=True)
    httpd = serve()
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch()
            page = browser.new_page(viewport={"width": 390, "height": 844},
                                    device_scale_factor=2, is_mobile=True,
                                    has_touch=True)
            for rel, name in SHOTS:
                page.goto(f"http://127.0.0.1:{PORT}/{rel}", wait_until="load")
                page.wait_for_timeout(400)
                page.screenshot(path=os.path.join(OUT, f"{name}.png"), full_page=False)
                print("shot", name)
            browser.close()
    finally:
        httpd.shutdown()


if __name__ == "__main__":
    main()
