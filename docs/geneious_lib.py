"""
geneious_lib.py -- read/write Geneious documents and place primer annotations.

A `.geneious` file is a ZIP archive containing a single UTF-8 XML entry.
The authoritative sequence lives in <charSequence>; annotations live inside
<sequenceAnnotations> as <annotation> blocks. Coordinates (minimumIndex /
maximumIndex) are 1-based, inclusive.

This module is deliberately dependency-free (stdlib only) so it can be reused
anywhere. See geneious_annot.py for the command-line front end.
"""

import csv
import io
import re
import zipfile
from xml.sax.saxutils import escape

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    """Reverse complement of a DNA string."""
    return seq.translate(COMPLEMENT)[::-1]


# --------------------------------------------------------------------------- #
# Container I/O  (.geneious  <->  inner XML)
# --------------------------------------------------------------------------- #
def read_geneious(path):
    """Return (inner_entry_name, xml_text) from a .geneious zip."""
    with zipfile.ZipFile(path) as z:
        name = z.namelist()[0]
        xml = z.read(name).decode("utf-8")
    return name, xml


def write_geneious(path, inner_name, xml_text):
    """Write xml_text into a new .geneious zip at `path` under `inner_name`."""
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(inner_name, xml_text.encode("utf-8"))


# --------------------------------------------------------------------------- #
# Parsing helpers
# --------------------------------------------------------------------------- #
def get_charsequence(xml_text):
    """Return the authoritative residue string from <charSequence>."""
    m = re.search(r"<charSequence>([A-Za-z]+)</charSequence>", xml_text)
    if not m:
        raise ValueError("no <charSequence> found")
    return m.group(1).upper()


def iter_annotations(xml_text):
    """Yield dicts summarising each <annotation> block (read-only view)."""
    for blk in re.findall(r"<annotation>.*?</annotation>", xml_text, re.S):
        desc = re.search(r"<description>(.*?)</description>", blk, re.S)
        typ = re.search(r"<type>(.*?)</type>", blk, re.S)
        iv = re.search(
            r"<minimumIndex>(\d+)</minimumIndex>"
            r"<maximumIndex>(\d+)</maximumIndex>"
            r"<direction>(\w+)</direction>",
            blk,
        )
        quals = {
            n: v
            for n, v in re.findall(
                r"<qualifier><name>(.*?)</name><value>(.*?)</value></qualifier>",
                blk,
                re.S,
            )
        }
        yield {
            "name": desc.group(1) if desc else "",
            "type": typ.group(1) if typ else "",
            "min": int(iv.group(1)) if iv else None,
            "max": int(iv.group(2)) if iv else None,
            "direction": iv.group(3) if iv else None,
            "qualifiers": quals,
        }


# --------------------------------------------------------------------------- #
# Primer alignment
# --------------------------------------------------------------------------- #
def find_binding_site(template, primer, min_anchor=8, strand=None):
    """Locate where `primer` anneals to `template`, anchored at the primer 3' end.

    Returns a dict, or None if no perfect anneal of length >= min_anchor exists.

    The primer's 3' end must match perfectly (that's what a polymerase needs);
    any non-matching bases are taken to be a 5' overhang/extension.

      strand 'forward': primer == top strand; 5' overhang is the left of S.
      strand 'reverse': primer == revcomp(top strand); 5' overhang is right of S.

    Result keys:
      strand       'forward' | 'reverse'
      start, end   1-based inclusive interval on the template (binding region)
      direction    'leftToRight' | 'rightToLeft'  (Geneious convention)
      binding_seq  the 3' portion of the primer that anneals (5'->3')
      extension    the 5' overhang of the primer that does NOT anneal (5'->3')
      occurrences  how many times the binding fragment occurs in the template
    """
    template = template.upper()
    primer = primer.upper()
    strands = ["forward", "reverse"] if strand is None else [strand]
    candidates = []

    for st in strands:
        S = primer if st == "forward" else revcomp(primer)
        n = len(S)
        # Scan from the longest possible anneal down to min_anchor, peeling
        # bases off the 5' overhang end until the remainder is found verbatim.
        for trim in range(0, n - min_anchor + 1):
            frag = S[trim:] if st == "forward" else S[: n - trim]
            pos = template.find(frag)
            if pos != -1:
                candidates.append(
                    {
                        "strand": st,
                        "pos": pos,
                        "blen": len(frag),
                        "overhang_len": trim,
                        "count": template.count(frag),
                    }
                )
                break  # first hit == longest anneal for this strand

    if not candidates:
        return None

    # Prefer the longest perfect anneal; tie-break toward forward strand.
    candidates.sort(key=lambda c: (c["blen"], c["strand"] == "forward"), reverse=True)
    best = candidates[0]
    overlen = best["overhang_len"]
    return {
        "strand": best["strand"],
        "start": best["pos"] + 1,
        "end": best["pos"] + best["blen"],
        "direction": "leftToRight" if best["strand"] == "forward" else "rightToLeft",
        "binding_seq": primer[overlen:],   # 3' annealing part
        "extension": primer[:overlen],     # 5' overhang (may be "")
        "occurrences": best["count"],
    }


def find_all_binding_sites(template, primer, min_anchor=15, strand=None):
    """Find EVERY locus where `primer` anneals with a perfect 3'-anchored run.

    A site is reported only if the contiguous, perfectly-matching stretch that
    INCLUDES the primer's 3'-terminal base is >= min_anchor bp. This is the
    guard against false positives: a primer whose 3' end matches only a few
    bases is rejected outright rather than dumping most of itself into an
    'extension'. Any non-matching 5' remainder becomes the 5' overhang.

    Returns a list of site dicts (possibly empty, or >1 for repeats), each:
      strand, start, end (1-based incl), direction, binding_seq, extension.
    The primer's 3'-terminal base always lies inside [start, end].
    """
    template = template.upper()
    primer = primer.upper()
    n = len(primer)
    strands = ["forward", "reverse"] if strand is None else [strand]
    sites = []

    for st in strands:
        S = primer if st == "forward" else revcomp(primer)
        m = len(S)
        if st == "forward":
            # primer 3' base == S[-1]; binding extends leftward (lower coords).
            anchor = S[-1]
            for e in range(len(template)):
                if template[e] != anchor:
                    continue
                L = 0
                while L < m and e - L >= 0 and template[e - L] == S[m - 1 - L]:
                    L += 1
                if L >= min_anchor:
                    sites.append({
                        "strand": "forward",
                        "start": e - L + 2,         # 1-based inclusive
                        "end": e + 1,
                        "direction": "leftToRight",
                        "binding_seq": primer[n - L:],
                        "extension": primer[:n - L],
                    })
        else:
            # primer 3' base == S[0]; binding extends rightward (higher coords).
            anchor = S[0]
            tn = len(template)
            for s in range(tn):
                if template[s] != anchor:
                    continue
                L = 0
                while L < m and s + L < tn and template[s + L] == S[L]:
                    L += 1
                if L >= min_anchor:
                    sites.append({
                        "strand": "reverse",
                        "start": s + 1,             # 1-based inclusive
                        "end": s + L,
                        "direction": "rightToLeft",
                        "binding_seq": primer[n - L:],
                        "extension": primer[:n - L],
                    })

    sites.sort(key=lambda x: (x["start"], x["strand"]))
    return sites


def parse_primers(text):
    """Parse (name, sequence) pairs from CSV/TSV text.

    Delimiter is auto-sniffed (tab vs comma). Column order is auto-detected:
    the cell that is pure DNA (ACGTN) is the sequence, the other is the name.
    Rows without a DNA cell (e.g. a header) are skipped.
    """
    delim = "\t" if text.count("\t") >= text.count(",") else ","
    rows = list(csv.reader(text.splitlines(), delimiter=delim))

    primers = []
    dna = re.compile(r"^[ACGTNacgtn]+$")
    for row in rows:
        cells = [c.strip() for c in row if c.strip()]
        if len(cells) < 2:
            continue
        seq = next((c for c in cells if dna.match(c)), None)
        if seq is None:
            continue  # header / comment / junk row
        name = next((c for c in cells if c != seq), cells[0])
        primers.append((name, seq.upper()))
    return primers


def read_primer_table(path):
    """Read (name, sequence) pairs from a CSV/TSV file. See parse_primers()."""
    with open(path, newline="") as f:
        return parse_primers(f.read())


# --------------------------------------------------------------------------- #
# Generating a fresh .geneious document from a raw sequence
# --------------------------------------------------------------------------- #
IUPAC = "ACGTUNRYSWKMBDHV"  # accepted residue alphabet (upper-cased)


def clean_sequence(text):
    """Extract a bare residue string from raw text or FASTA.

    Drops FASTA header lines (starting with '>'), whitespace and any character
    outside the IUPAC nucleotide alphabet. Returns an upper-case string.
    """
    out = []
    for line in text.splitlines():
        if line.startswith(">"):
            continue
        out.append(line)
    seq = "".join(out).upper()
    return "".join(c for c in seq if c in IUPAC)


def make_geneious(name, sequence, topology="linear", annotations_xml=""):
    """Build a minimal but valid Geneious nucleotide document XML string.

    Contains the load-bearing structure only (no database-linkage, history, or
    view-state cruft): a GenBankNucleotideSequence whose fields and
    originalElement carry the sequence, name and (optional) annotations.
    """
    seq = clean_sequence(sequence)
    if not seq:
        raise ValueError("no valid nucleotide residues found in input")
    L = len(seq)
    gc = (seq.count("G") + seq.count("C")) / L * 100.0
    nm = escape(name or "sequence")
    ts = "1700000000000"  # fixed placeholder epoch-millis; Geneious resets on save

    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<geneious version="2026.1.2">\n'
        '<versions version="2026.1.2" minimumVersion="7.1"/>\n'
        '<geneiousDocument '
        'class="com.biomatters.plugins.ncbi.documents.GenBankNucleotideSequence" '
        'version="1.3-11" revisionNumber="1" geneiousVersion="2026.1.2" '
        'geneiousVersionMinimum="7.1" PluginDocument_FormatLastChanged="7.1" '
        'PluginDocument_FormatLastExtended="7.1" '
        'PluginDocument_OldestVersionSerializableTo="6.0" isReferenceOnly="false">'
        "<hiddenFields>"
        "<description></description>"
        "<cache_name>%s</cache_name>"
        '<override_modified_date type="date">%s</override_modified_date>'
        '<cache_created type="date">%s</cache_created>'
        "</hiddenFields>"
        "<fields>"
        "<molType>DNA</molType>"
        "<topology>%s</topology>"
        '<sequence_length type="int">%d</sequence_length>'
        "<sequence_residues>%s</sequence_residues>"
        "<geneticCode>Standard</geneticCode>"
        '<gcPercent decimalPlacesDisplayed="1" type="percent">%s</gcPercent>'
        '<modified_date type="date">%s</modified_date>'
        "</fields>"
        "<originalElement><XMLSerialisableRootElement>"
        "<fields>"
        "<geneticCode>Standard</geneticCode>"
        "<topology>%s</topology>"
        "<molType>DNA</molType>"
        "</fields>"
        '<urn type="urn">urn:sequence:local:generated</urn>'
        '<created type="date">%s</created>'
        "<storedFields>"
        '<standardField code="molType" />'
        '<standardField code="topology" />'
        '<standardField code="geneticCode" />'
        "</storedFields>"
        "<name>%s</name>"
        "<description />"
        "<sequenceAnnotations>%s</sequenceAnnotations>"
        "<charSequence>%s</charSequence>"
        "<INSD_originalElements>"
        "<INSDSeq_length>%d</INSDSeq_length>"
        "<INSDSeq_locus>%s</INSDSeq_locus>"
        "<INSDSeq_strandedness>double</INSDSeq_strandedness>"
        "</INSD_originalElements>"
        "</XMLSerialisableRootElement></originalElement>"
        "</geneiousDocument>\n"
        "</geneious>\n"
    ) % (
        nm, ts, ts,
        escape(topology), L, seq, ("%.1f" % gc), ts,
        escape(topology), ts, nm, annotations_xml, seq, L, nm,
    )


# --------------------------------------------------------------------------- #
# Annotation construction / insertion
# --------------------------------------------------------------------------- #
def build_primer_annotation(
    name, start, end, direction, binding_seq, extension="",
    mismatches=0, created_by="manual", extra_qualifiers=None,
):
    """Return an <annotation> XML string for a primer_bind feature.

    Mirrors the Geneious '5oh' representation: the annotation spans only the
    annealing region; a 5' overhang is stored in the 'Extension' qualifier.
    Full oligo (5'->3') == Extension + Sequence.
    """
    quals = [("Sequence", binding_seq)]
    if extension:
        quals.append(("Extension", extension))
    quals.append(("Mismatches", str(mismatches)))
    if created_by:
        quals.append(("created by", created_by))
    if extra_qualifiers:
        quals.extend(extra_qualifiers)

    qx = "".join(
        "<qualifier><name>%s</name><value>%s</value></qualifier>"
        % (escape(n), escape(v))
        for n, v in quals
    )
    return (
        "<annotation>"
        "<description>%s</description>"
        "<type>primer_bind</type>"
        "<intervals><interval>"
        "<minimumIndex>%d</minimumIndex>"
        "<maximumIndex>%d</maximumIndex>"
        "<direction>%s</direction>"
        "</interval></intervals>"
        "<qualifiers>%s</qualifiers>"
        "</annotation>"
    ) % (escape(name), start, end, direction, qx)


def insert_annotation(xml_text, annotation_xml):
    """Insert an <annotation> block just before </sequenceAnnotations>."""
    marker = "</sequenceAnnotations>"
    if xml_text.count(marker) != 1:
        raise ValueError("expected exactly one </sequenceAnnotations>, found %d"
                         % xml_text.count(marker))
    return xml_text.replace(marker, annotation_xml + marker)


# --------------------------------------------------------------------------- #
# High-level: two inputs (sequence + primers) -> annotated .geneious bytes
# Shared by the local server (geneious_app.py) and the in-browser app (docs/).
# --------------------------------------------------------------------------- #
def looks_like_geneious(data):
    """(is_geneious, inner_xml) for a bytes blob that may be a .geneious zip."""
    try:
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            inner = z.read(z.namelist()[0]).decode("utf-8", "replace")
        return ("<geneious" in inner and "<charSequence>" in inner), inner
    except Exception:
        return False, ""


def safe_filename(name):
    """Sanitise a document name into a '<name>.geneious' filename."""
    keep = "-_." + "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    base = "".join(c if c in keep else "_" for c in (name or "output")).strip("_")
    return (base or "output") + ".geneious"


def fasta_name(text):
    """First FASTA header token, or None if the text has no '>' header."""
    for line in text.splitlines():
        if line.startswith(">"):
            head = line[1:].strip()
            return head.split()[0] if head else None
    return None


def build_document(seq_source, primers_text, seq_name="", min_anchor=15):
    """Turn a sequence source + primer table into an annotated .geneious.

    seq_source is one of:
        {"kind": "paste", "text": <raw sequence or FASTA>}
        {"kind": "file",  "bytes": <bytes of a .geneious OR FASTA/text file>}

    A .geneious blob is annotated in place (existing annotations kept); raw
    text / FASTA generates a fresh document. Returns (filename, bytes, log).
    """
    if seq_source.get("kind") == "file":
        data = seq_source["bytes"]
        is_gen, inner = looks_like_geneious(data)
        if is_gen:
            base_xml = inner
            template = get_charsequence(base_xml)
            doc_name = seq_name or "annotated"
        else:
            text = data.decode("utf-8", "replace")
            doc_name = seq_name or fasta_name(text) or "sequence"
            template = clean_sequence(text)
            base_xml = make_geneious(doc_name, text)
    else:
        text = seq_source.get("text", "")
        doc_name = seq_name or fasta_name(text) or "sequence"
        template = clean_sequence(text)
        base_xml = make_geneious(doc_name, text)

    if not template:
        raise ValueError("the target sequence is empty / contains no valid residues")

    primers = parse_primers(primers_text)
    log = ["Target: %s  (%d bp)" % (doc_name, len(template)),
           "Primers loaded: %d   |   anchor cutoff: >= %d bp\n" % (len(primers), min_anchor)]

    xml = base_xml
    matched = multi = 0
    for name, seq in primers:
        sites = find_all_binding_sites(template, seq, min_anchor=min_anchor)
        if not sites:
            log.append("  %-16s NO MATCH (>=%d bp) -- skipped" % (name, min_anchor))
            continue
        matched += 1
        is_multi = len(sites) > 1
        if is_multi:
            multi += 1
            log.append("  %-16s *** MULTIPLE (%d) matches -- annotating all ***"
                       % (name, len(sites)))
        for i, s in enumerate(sites, 1):
            ann_name = "%s_%d" % (name, i) if is_multi else name
            ext = s["extension"]
            xml = insert_annotation(xml, build_primer_annotation(
                ann_name, s["start"], s["end"], s["direction"],
                s["binding_seq"], ext, mismatches=0, created_by="geneious_tools"))
            log.append("  %-16s %-7s %d-%d  bind=%dnt%s" % (
                ann_name, s["strand"], s["start"], s["end"], len(s["binding_seq"]),
                ("  ext(5')=%s" % ext) if ext else "  (full match, no ext)"))

    log.append("\nDone: %d/%d primers matched%s." % (
        matched, len(primers), (", %d multi-locus" % multi) if multi else ""))

    filename = safe_filename(doc_name)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(filename, xml.encode("utf-8"))
    return filename, buf.getvalue(), "\n".join(log)
