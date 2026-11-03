"""Check the rendered site at phone widths, using a real headless browser.

Loads each page at 390 px and 360 px and asserts nothing overflows the viewport
horizontally. Prints the widest offending element when it does.
"""
import json
import os
import sys
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

ROOT = r"D:\site\docs"
PAGES = ["index.html", "about.html", "articles/outsides.html", "articles/models.html",
         "articles/signal.html", "articles/systems.html"]
WIDTHS = [390, 360]
PORT = 8899

PROBE = """
(function () {
  var out = {scrollW: document.documentElement.scrollWidth, vw: window.innerWidth, overflow: []};
  document.querySelectorAll('body *').forEach(function (n) {
    var r = n.getBoundingClientRect();
    if (r.width > 0 && r.right > out.vw + 1) {
      out.overflow.push({el: n.tagName + '.' + (n.className || ''), right: Math.round(r.right),
                         w: Math.round(r.width)});
    }
  });
  out.overflow = out.overflow.slice(0, 6);
  return out;
})()
"""


def serve():
    handler = partial(SimpleHTTPRequestHandler, directory=ROOT)
    httpd = ThreadingHTTPServer(("127.0.0.1", PORT), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    return httpd


def main():
    from playwright.sync_api import sync_playwright

    httpd = serve()
    failures = 0
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch()
            for width in WIDTHS:
                page = browser.new_page(viewport={"width": width, "height": 800})
                for rel in PAGES:
                    page.goto(f"http://127.0.0.1:{PORT}/{rel}", wait_until="load")
                    page.wait_for_timeout(250)
                    res = page.evaluate(PROBE)
                    bad = res["scrollW"] > width + 1 or res["overflow"]
                    status = "FAIL" if bad else "ok"
                    if bad:
                        failures += 1
                    print(f"{width:4}px  {rel:28} {status:5} scrollW={res['scrollW']}")
                    for o in res["overflow"]:
                        print(f"          overflow: {o['el'][:56]} right={o['right']} w={o['w']}")
                page.close()
            browser.close()
    finally:
        httpd.shutdown()

    print(f"\n{failures} failing page/width combinations")
    assert failures == 0, "mobile layout overflows"


if __name__ == "__main__":
    sys.exit(main())
