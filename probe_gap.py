"""Find which element is eating the vertical space above the page content on a phone."""
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

from playwright.sync_api import sync_playwright

ROOT = r"D:\site\docs"
PORT = 8897

PROBE = """
(function () {
  var out = [];
  var nav = document.querySelector('.navbar');
  var out_r = nav ? nav.getBoundingClientRect() : null;
  var walk = document.querySelector('#quarto-document-content') || document.body;
  var first = walk.querySelector('p, .quarto-title-block, .panel-tabset');
  var r = first ? first.getBoundingClientRect() : null;
  var chain = [];
  var el = first;
  while (el && el !== document.body) {
    var cs = getComputedStyle(el);
    chain.push({el: el.tagName + '.' + (el.className || '').toString().slice(0,40),
                mt: cs.marginTop, pt: cs.paddingTop, top: Math.round(el.getBoundingClientRect().top)});
    el = el.parentElement;
  }
  return {navBottom: out_r ? Math.round(out_r.bottom) : null,
          firstTop: r ? Math.round(r.top) : null, chain: chain};
})()
"""

handler = partial(SimpleHTTPRequestHandler, directory=ROOT)
httpd = ThreadingHTTPServer(("127.0.0.1", PORT), handler)
threading.Thread(target=httpd.serve_forever, daemon=True).start()

with sync_playwright() as p:
    b = p.chromium.launch()
    pg = b.new_page(viewport={"width": 390, "height": 844})
    for rel in ["articles/systems.html", "index.html"]:
        pg.goto(f"http://127.0.0.1:{PORT}/{rel}", wait_until="load")
        pg.wait_for_timeout(300)
        res = pg.evaluate(PROBE)
        print(f"\n--- {rel} --- navBottom={res['navBottom']} firstTop={res['firstTop']}")
        for c in res["chain"]:
            print(f"   {c['el'][:44]:46} top={c['top']:4}  margin-top={c['mt']:10} padding-top={c['pt']}")
    b.close()
httpd.shutdown()
