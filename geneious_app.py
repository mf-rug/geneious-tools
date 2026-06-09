#!/usr/bin/env python3
"""
geneious_app.py -- minimal local web GUI for the Geneious primer annotator.

Run it, open the browser tab it prints, paste (or pick a file for) a target
sequence and a primer table, and download the annotated .geneious file.

    python3 geneious_app.py            # serves http://127.0.0.1:8765
    python3 geneious_app.py 9000       # custom port

Dependency-free (Python standard library only). Files are read in the browser
and sent as base64, so the server never touches the filesystem.
"""

import base64
import io
import json
import sys
import webbrowser
import zipfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import geneious_lib as gl


# --------------------------------------------------------------------------- #
# Core: turn the two inputs into an annotated .geneious (bytes) + a log
# --------------------------------------------------------------------------- #
def _looks_like_geneious(data):
    try:
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            inner = z.read(z.namelist()[0]).decode("utf-8", "replace")
        return "<geneious" in inner and "<charSequence>" in inner, inner
    except Exception:
        return False, ""


def _safe_filename(name):
    keep = "-_." + "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    base = "".join(c if c in keep else "_" for c in (name or "output")).strip("_")
    return (base or "output") + ".geneious"


def build_document(seq_source, primers_text, seq_name, min_anchor):
    """seq_source is either {'kind':'paste','text':..} or {'kind':'file','bytes':..}.

    Returns (filename, geneious_bytes, log_text).
    """
    # 1. Resolve the base document + template sequence.
    if seq_source["kind"] == "file":
        data = seq_source["bytes"]
        is_gen, inner = _looks_like_geneious(data)
        if is_gen:
            base_xml = inner
            template = gl.get_charsequence(base_xml)
            doc_name = seq_name or "annotated"
        else:  # treat the uploaded file as raw text / FASTA
            text = data.decode("utf-8", "replace")
            doc_name = seq_name or _fasta_name(text) or "sequence"
            template = gl.clean_sequence(text)
            base_xml = gl.make_geneious(doc_name, template)
    else:  # pasted raw sequence
        text = seq_source["text"]
        doc_name = seq_name or _fasta_name(text) or "sequence"
        template = gl.clean_sequence(text)
        base_xml = gl.make_geneious(doc_name, template)

    if not template:
        raise ValueError("the target sequence is empty / contains no valid residues")

    # 2. Parse primers and annotate.
    primers = gl.parse_primers(primers_text)
    log = []
    log.append("Target: %s  (%d bp)" % (doc_name, len(template)))
    log.append("Primers loaded: %d   |   anchor cutoff: >= %d bp\n" % (len(primers), min_anchor))

    xml = base_xml
    matched = multi = 0
    for name, seq in primers:
        sites = gl.find_all_binding_sites(template, seq, min_anchor=min_anchor)
        if not sites:
            log.append("  %-16s NO MATCH (>=%d bp) -- skipped" % (name, min_anchor))
            continue
        matched += 1
        is_multi = len(sites) > 1
        if is_multi:
            multi += 1
            log.append("  %-16s *** MULTIPLE (%d) matches -- annotating all ***" % (name, len(sites)))
        for i, s in enumerate(sites, 1):
            ann_name = "%s_%d" % (name, i) if is_multi else name
            ext = s["extension"]
            xml = gl.insert_annotation(
                xml,
                gl.build_primer_annotation(
                    ann_name, s["start"], s["end"], s["direction"],
                    s["binding_seq"], ext, mismatches=0, created_by="geneious_app",
                ),
            )
            log.append("  %-16s %-7s %d-%d  bind=%dnt%s" % (
                ann_name, s["strand"], s["start"], s["end"], len(s["binding_seq"]),
                ("  ext(5')=%s" % ext) if ext else "  (full match, no ext)"))

    log.append("\nDone: %d/%d primers matched%s." % (
        matched, len(primers), (", %d multi-locus" % multi) if multi else ""))

    # 3. Zip it up.
    filename = _safe_filename(doc_name)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(filename, xml.encode("utf-8"))
    return filename, buf.getvalue(), "\n".join(log)


def _fasta_name(text):
    for line in text.splitlines():
        if line.startswith(">"):
            return line[1:].strip().split()[0] if line[1:].strip() else None
    return None


# --------------------------------------------------------------------------- #
# HTTP handler
# --------------------------------------------------------------------------- #
class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):  # quiet console
        pass

    def do_GET(self):
        if self.path not in ("/", "/index.html"):
            self.send_error(404)
            return
        body = PAGE.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        if self.path != "/generate":
            self.send_error(404)
            return
        length = int(self.headers.get("Content-Length", 0))
        try:
            req = json.loads(self.rfile.read(length).decode("utf-8"))
            seq_source = _decode_source(req, "seq")
            primers_text = _decode_text(req, "primer")
            if not primers_text.strip():
                raise ValueError("no primers provided")
            fname, data, log = build_document(
                seq_source, primers_text,
                req.get("seq_name", "").strip(),
                int(req.get("min_anchor", 15)),
            )
            payload = {
                "ok": True,
                "filename": fname,
                "file_b64": base64.b64encode(data).decode("ascii"),
                "log": log,
            }
        except Exception as e:
            payload = {"ok": False, "error": "%s: %s" % (type(e).__name__, e)}

        out = json.dumps(payload).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(out)))
        self.end_headers()
        self.wfile.write(out)


def _decode_source(req, prefix):
    if req.get(prefix + "_mode") == "file":
        b64 = req.get(prefix + "_file_b64", "")
        if not b64:
            raise ValueError("no sequence file selected")
        return {"kind": "file", "bytes": base64.b64decode(b64)}
    return {"kind": "paste", "text": req.get(prefix + "_text", "")}


def _decode_text(req, prefix):
    if req.get(prefix + "_mode") == "file":
        b64 = req.get(prefix + "_file_b64", "")
        return base64.b64decode(b64).decode("utf-8", "replace") if b64 else ""
    return req.get(prefix + "_text", "")


# --------------------------------------------------------------------------- #
# The single-page UI
# --------------------------------------------------------------------------- #
PAGE = r"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Geneious primer annotator</title>
<style>
  :root { --bg:#0f1419; --panel:#1a2129; --line:#2c3744; --fg:#e6edf3; --mut:#8b98a5;
          --accent:#3fb950; --accent2:#2ea043; }
  * { box-sizing:border-box; }
  body { margin:0; background:var(--bg); color:var(--fg);
         font:14px/1.5 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; }
  .wrap { max-width:860px; margin:0 auto; padding:28px 20px 60px; }
  h1 { font-size:20px; margin:0 0 4px; }
  .sub { color:var(--mut); margin:0 0 24px; }
  .card { background:var(--panel); border:1px solid var(--line); border-radius:10px;
          padding:18px; margin-bottom:18px; }
  .card h2 { font-size:14px; margin:0 0 12px; letter-spacing:.04em; text-transform:uppercase;
             color:var(--mut); }
  .tabs { display:flex; gap:6px; margin-bottom:10px; }
  .tab { padding:5px 12px; border:1px solid var(--line); border-radius:6px; cursor:pointer;
         color:var(--mut); background:transparent; font:inherit; }
  .tab.on { color:var(--fg); border-color:var(--accent); background:#13311c; }
  textarea { width:100%; min-height:120px; resize:vertical; background:#0d1117; color:var(--fg);
             border:1px solid var(--line); border-radius:6px; padding:10px; font:inherit; }
  input[type=text], input[type=number] { background:#0d1117; color:var(--fg);
             border:1px solid var(--line); border-radius:6px; padding:8px 10px; font:inherit; }
  input[type=file]{ color:var(--mut); font:inherit; }
  .row { display:flex; gap:18px; flex-wrap:wrap; align-items:flex-end; }
  .row > div { flex:1; min-width:180px; }
  label.fld { display:block; color:var(--mut); margin-bottom:5px; font-size:12px; }
  .hide { display:none; }
  button.go { background:var(--accent2); color:#fff; border:0; border-radius:8px;
              padding:12px 22px; font:inherit; font-weight:600; cursor:pointer; }
  button.go:hover { background:var(--accent); }
  button.go:disabled { opacity:.5; cursor:default; }
  pre { background:#0d1117; border:1px solid var(--line); border-radius:8px; padding:14px;
        white-space:pre-wrap; word-break:break-word; max-height:340px; overflow:auto; }
  .err { color:#ff7b72; }
  a.dl { color:var(--accent); }
  .hint { color:var(--mut); font-size:12px; margin-top:6px; }
</style></head>
<body><div class="wrap">
  <h1>Geneious primer annotator</h1>
  <p class="sub">Paste or upload a target sequence + primer table → download an annotated .geneious file.</p>

  <div class="card">
    <h2>1 · Target sequence</h2>
    <div class="tabs" data-group="seq">
      <button class="tab on" data-mode="paste">Paste</button>
      <button class="tab" data-mode="file">File</button>
    </div>
    <div data-pane="seq-paste">
      <textarea id="seq_text" placeholder="Raw sequence or FASTA, e.g.&#10;>my_plasmid&#10;ATGCAT...&#10;(boilerplate for a valid .geneious is generated automatically)"></textarea>
    </div>
    <div data-pane="seq-file" class="hide">
      <input type="file" id="seq_file" accept=".geneious,.fasta,.fa,.fna,.txt,.seq">
      <div class="hint">A .geneious file is annotated in place (existing annotations kept). A FASTA/text file becomes a new .geneious.</div>
    </div>
    <div class="row" style="margin-top:12px">
      <div><label class="fld">Sequence name (optional)</label><input type="text" id="seq_name" placeholder="sequence"></div>
      <div style="flex:0 0 160px"><label class="fld">Min 3' anneal (bp)</label><input type="number" id="min_anchor" value="15" min="1"></div>
    </div>
  </div>

  <div class="card">
    <h2>2 · Primers (CSV / TSV: name, sequence)</h2>
    <div class="tabs" data-group="primer">
      <button class="tab on" data-mode="paste">Paste</button>
      <button class="tab" data-mode="file">File</button>
    </div>
    <div data-pane="primer-paste">
      <textarea id="primer_text" placeholder="name,sequence&#10;fwd1,AGGTCCCCGAAGCTGCTATTTCACG&#10;rev1,AATGAATGGTTAGCCCATCATCTCTTC"></textarea>
    </div>
    <div data-pane="primer-file" class="hide">
      <input type="file" id="primer_file" accept=".csv,.tsv,.txt">
    </div>
  </div>

  <button class="go" id="go">Generate .geneious</button>
  <div id="out" style="margin-top:20px"></div>
</div>

<script>
const modes = {seq:"paste", primer:"paste"};
document.querySelectorAll(".tabs").forEach(t => {
  const grp = t.dataset.group;
  t.querySelectorAll(".tab").forEach(b => b.onclick = () => {
    modes[grp] = b.dataset.mode;
    t.querySelectorAll(".tab").forEach(x => x.classList.toggle("on", x === b));
    document.querySelector(`[data-pane="${grp}-paste"]`).classList.toggle("hide", b.dataset.mode!=="paste");
    document.querySelector(`[data-pane="${grp}-file"]`).classList.toggle("hide", b.dataset.mode!=="file");
  });
});

function readB64(input) {
  return new Promise((res, rej) => {
    const f = input.files[0];
    if (!f) return res(null);
    const r = new FileReader();
    r.onload = () => res(r.result.split(",")[1]);   // strip data: prefix
    r.onerror = rej;
    r.readAsDataURL(f);
  });
}

document.getElementById("go").onclick = async () => {
  const btn = document.getElementById("go"), out = document.getElementById("out");
  btn.disabled = true; out.innerHTML = "<pre>Working…</pre>";
  try {
    const req = {
      seq_mode: modes.seq, primer_mode: modes.primer,
      seq_text: document.getElementById("seq_text").value,
      primer_text: document.getElementById("primer_text").value,
      seq_name: document.getElementById("seq_name").value,
      min_anchor: document.getElementById("min_anchor").value || 15,
      seq_file_b64: modes.seq==="file" ? await readB64(document.getElementById("seq_file")) : null,
      primer_file_b64: modes.primer==="file" ? await readB64(document.getElementById("primer_file")) : null,
    };
    const resp = await fetch("/generate", {method:"POST", headers:{"Content-Type":"application/json"}, body:JSON.stringify(req)});
    const j = await resp.json();
    if (!j.ok) { out.innerHTML = `<pre class="err">${j.error}</pre>`; return; }
    const bytes = Uint8Array.from(atob(j.file_b64), c => c.charCodeAt(0));
    const url = URL.createObjectURL(new Blob([bytes], {type:"application/octet-stream"}));
    out.innerHTML = `<p><a class="dl" id="dl" href="${url}" download="${j.filename}">⤓ Download ${j.filename}</a></p><pre>${j.log.replace(/</g,"&lt;")}</pre>`;
    document.getElementById("dl").click();   // auto-trigger download
  } catch (e) {
    out.innerHTML = `<pre class="err">${e}</pre>`;
  } finally { btn.disabled = false; }
};
</script>
</body></html>
"""


def main():
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8765
    url = "http://127.0.0.1:%d" % port
    srv = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    print("Geneious primer annotator -> %s   (Ctrl-C to stop)" % url)
    try:
        webbrowser.open(url)
    except Exception:
        pass
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped.")


if __name__ == "__main__":
    main()
