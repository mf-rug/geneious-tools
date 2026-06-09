# geneious-tools

Read [Geneious](https://www.geneious.com/) documents and programmatically add
primer annotations — without opening Geneious.

A `.geneious` file is just a ZIP wrapping a single UTF-8 XML document. The
sequence lives in `<charSequence>`; annotations are `<annotation>` blocks inside
`<sequenceAnnotations>`, with 1-based inclusive coordinates. These tools parse
that format, locate where primers anneal, and write back valid `.geneious`
files that Geneious opens normally.

Pure Python standard library — no dependencies.

## Features

- **Batch primer annotation** from a CSV/TSV of `name, sequence`.
- **Both strands** searched automatically, anchored at the primer's 3' end.
- **3' match + 5' overhang handling.** A primer that matches fully is annotated
  with no extension; a primer whose 3' end matches but whose 5' end does not has
  everything from the first mismatch stored as a 5' `Extension` qualifier
  (full oligo = `Extension` + `Sequence`).
- **Strict anneal cutoff** (default ≥15 bp of perfect, contiguous 3' match).
  Primers whose 3' end matches only a few bases are rejected rather than dumping
  most of the oligo into a bogus "extension".
- **Repeat-aware.** If a primer matches multiple loci, every match is annotated
  (`name_1`, `name_2`, …) and the multi-match is flagged.
- **Non-destructive.** The input file is never modified; results go to a new file.

## Web app (nothing to install)

**→ https://mf-rug.github.io/geneious-tools/**

Runs entirely in your browser via [Pyodide](https://pyodide.org) (Python compiled
to WebAssembly). No install, no account, no server — and your sequences never
leave your machine. First load fetches the ~7 MB runtime (cached afterwards).

It gives you two inputs — **target sequence** and **primers** — each accepting
either a **pasted value** or a **file**:

- *Target sequence* — paste raw sequence/FASTA (valid `.geneious` boilerplate is
  generated automatically), or upload a file. Uploading an existing `.geneious`
  annotates it in place (keeping its annotations); uploading FASTA/text makes a
  new document.
- *Primers* — paste or upload a CSV/TSV of `name, sequence`.

Click **Generate** and the annotated `.geneious` downloads, with a log of which
primers matched where.

## Local GUI (same app, offline)

The same interface also runs as a tiny local server (no internet needed):

```bash
python3 geneious_app.py            # opens http://127.0.0.1:8765
python3 geneious_app.py 9000       # custom port
```

It gives you two inputs — **target sequence** and **primers** — each accepting
either a **pasted value** or a **file**:

- *Target sequence* — paste raw sequence/FASTA (valid `.geneious` boilerplate is
  generated automatically), or upload a file. Uploading an existing `.geneious`
  annotates it in place (keeping its annotations); uploading FASTA/text makes a
  new document.
- *Primers* — paste or upload a CSV/TSV of `name, sequence`.

Click **Generate** and the annotated `.geneious` downloads. A log shows which
primers matched, where, and any multi-locus hits. Files are read in the browser
and sent as base64, so the server never touches your filesystem. Standard
library only — nothing to install.

## Command line

```bash
# Annotate a whole primer table onto a fresh copy
python3 geneious_annot.py batch  in.geneious  out.geneious  primers.tsv

# List existing annotations in a file
python3 geneious_annot.py list  in.geneious

# Add a single primer (auto strand + overhang detection)
python3 geneious_annot.py add-primer  in.geneious  out.geneious \
    --name myprimer --seq AATGAATGGTTAGCCCATCATCTCTTC
```

Useful flags for `batch` / `add-primer`:

- `--min-anchor N` — minimum perfect 3' anneal to accept a match (batch default 15).
- `--strand forward|reverse` — restrict to one strand (default: search both).
- `--created-by NAME` — value for the annotation's `created by` qualifier.

### Primer table

CSV or TSV with a `name` column and a `sequence` column. The delimiter and
column order are auto-detected, and a header row is auto-skipped. See
[`example_primers.tsv`](example_primers.tsv):

```
name	sequence
fwd_5oh	AGGTCCCCGAAGCTGCTATTTCACG
rev_5ext	AATGAATGGTTAGCCCATCATCTCTTC
```

## Library API

`geneious_lib.py` is usable directly:

```python
import geneious_lib as gl

inner, xml = gl.read_geneious("in.geneious")
template = gl.get_charsequence(xml)

for site in gl.find_all_binding_sites(template, "AATGAATGGTTAGCCCATCATCTCTTC", min_anchor=15):
    print(site)   # strand, start, end, direction, binding_seq, extension

ann = gl.build_primer_annotation("p1", site["start"], site["end"],
                                 site["direction"], site["binding_seq"], site["extension"])
gl.write_geneious("out.geneious", "out.geneious", gl.insert_annotation(xml, ann))
```

## Notes

- Coordinates (`minimumIndex` / `maximumIndex`) are 1-based and inclusive;
  `direction` is `leftToRight` (forward) or `rightToLeft` (reverse).
- The authoritative sequence is `<charSequence>`; some files also carry a
  truncated `<sequence_residues>` copy, which is ignored.
- Tested against Geneious documents from versions 8.x through 2026.x.

## Development

`geneious_lib.py` (repo root) is the single source of truth. The GitHub Pages
app serves its own copy at `docs/geneious_lib.py`; a versioned pre-commit hook
keeps the two byte-identical. After cloning, activate it once:

```bash
git config core.hooksPath .githooks
```

Then any commit that changes the root library auto-syncs and stages the `docs/`
copy — no manual copying, no drift.

## License

MIT — see [LICENSE](LICENSE).
