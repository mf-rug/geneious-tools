#!/usr/bin/env python3
"""
geneious_annot.py -- command-line tool to inspect Geneious documents and add
primer annotations (with automatic strand detection and 5' overhang handling).

Examples
--------
List existing annotations:
    python3 geneious_annot.py list  in.geneious

Add a primer (auto-detects strand; 3' must match, 5' overhang -> Extension):
    python3 geneious_annot.py add-primer in.geneious out.geneious \
        --name rev_ext --seq AATGAATGGTTAGCCCATCATCTCTTC

The input file is never modified; results are always written to a new file.
"""

import argparse
import sys

import geneious_lib as gl


def cmd_list(args):
    _, xml = gl.read_geneious(args.input)
    for a in gl.iter_annotations(xml):
        loc = (
            "%s..%s %s" % (a["min"], a["max"], a["direction"])
            if a["min"] is not None
            else "(no interval)"
        )
        seq = a["qualifiers"].get("Sequence", "")
        ext = a["qualifiers"].get("Extension", "")
        extra = (" +ext[%s]" % ext) if ext else ""
        print("%-22s %-12s %-22s %s%s" % (a["name"], a["type"], loc, seq, extra))


def cmd_add_primer(args):
    inner, xml = gl.read_geneious(args.input)
    template = gl.get_charsequence(xml)

    site = gl.find_binding_site(
        template, args.seq, min_anchor=args.min_anchor, strand=args.strand
    )
    if site is None:
        sys.exit(
            "ERROR: no perfect 3' anneal of length >= %d found on either strand."
            % args.min_anchor
        )

    # Report what was detected.
    print("Primer            : %s (%d nt)" % (args.seq, len(args.seq)))
    print("Strand            : %s (%s)" % (site["strand"], site["direction"]))
    print("Binding region    : %d-%d  (%d nt, %d mismatches)"
          % (site["start"], site["end"], len(site["binding_seq"]), 0))
    print("  binding (3')    : %s" % site["binding_seq"])
    print("  5' extension    : %s" % (site["extension"] or "(none)"))
    if site["occurrences"] > 1:
        print("  WARNING         : binding fragment occurs %d times; using the first."
              % site["occurrences"])

    # Sanity check: the annealing region's reverse/forward image must equal the
    # template substring at the chosen coordinates.
    region = template[site["start"] - 1 : site["end"]]
    expect = (
        site["binding_seq"]
        if site["strand"] == "forward"
        else gl.revcomp(site["binding_seq"])
    )
    assert region == expect, "internal coordinate mismatch: %s != %s" % (region, expect)

    ann = gl.build_primer_annotation(
        name=args.name,
        start=site["start"],
        end=site["end"],
        direction=site["direction"],
        binding_seq=site["binding_seq"],
        extension=site["extension"],
        mismatches=0,
        created_by=args.created_by,
    )
    new_xml = gl.insert_annotation(xml, ann)

    inner_out = args.inner_name or args.output.rsplit("/", 1)[-1]
    gl.write_geneious(args.output, inner_out, new_xml)

    # Verify the written copy round-trips and contains the new annotation.
    _, chk = gl.read_geneious(args.output)
    names = [a["name"] for a in gl.iter_annotations(chk)]
    ok = args.name in names
    print("Wrote             : %s  (annotations %d -> %d)"
          % (args.output, len(list(gl.iter_annotations(xml))), len(names)))
    print("Verification      : %s" % ("OK" if ok else "FAILED — annotation not found"))
    if not ok:
        sys.exit(1)


def cmd_batch(args):
    inner, xml = gl.read_geneious(args.input)
    template = gl.get_charsequence(xml)
    primers = gl.read_primer_table(args.primers)
    print("Template          : %d bp" % len(template))
    print("Primers loaded    : %d  (from %s)" % (len(primers), args.primers))
    print("Anchor cutoff     : >= %d bp perfect 3' anneal\n" % args.min_anchor)

    n_before = len(list(gl.iter_annotations(xml)))
    report = []  # (name, status, n_sites)

    for name, seq in primers:
        sites = gl.find_all_binding_sites(
            template, seq, min_anchor=args.min_anchor, strand=args.strand
        )
        if not sites:
            print("  %-16s NO MATCH (>=%d bp) -- skipped" % (name, args.min_anchor))
            report.append((name, "no_match", 0))
            continue

        multi = len(sites) > 1
        if multi:
            print("  %-16s *** MULTIPLE (%d) matches -- annotating all ***"
                  % (name, len(sites)))
        for i, site in enumerate(sites, 1):
            ann_name = "%s_%d" % (name, i) if multi else name
            ext = site["extension"]
            ann = gl.build_primer_annotation(
                name=ann_name,
                start=site["start"], end=site["end"],
                direction=site["direction"],
                binding_seq=site["binding_seq"], extension=ext,
                mismatches=0, created_by=args.created_by,
            )
            xml = gl.insert_annotation(xml, ann)
            print("  %-16s %-7s %d-%d  bind=%dnt%s"
                  % (ann_name, site["strand"], site["start"], site["end"],
                     len(site["binding_seq"]),
                     ("  ext(5')=%s" % ext) if ext else "  (full match, no ext)"))
        report.append((name, "multi" if multi else "ok", len(sites)))

    inner_out = args.inner_name or args.output.rsplit("/", 1)[-1]
    gl.write_geneious(args.output, inner_out, xml)
    n_after = len(list(gl.iter_annotations(xml)))

    print("\nSummary:")
    matched = sum(1 for _, s, _ in report if s != "no_match")
    multis = [n for n, s, _ in report if s == "multi"]
    nomatch = [n for n, s, _ in report if s == "no_match"]
    print("  primers matched : %d / %d" % (matched, len(report)))
    if multis:
        print("  MULTI-MATCH     : %s" % ", ".join(multis))
    if nomatch:
        print("  no match        : %s" % ", ".join(nomatch))
    print("  annotations     : %d -> %d  (+%d)" % (n_before, n_after, n_after - n_before))
    print("  wrote           : %s" % args.output)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    pl = sub.add_parser("list", help="list annotations in a .geneious file")
    pl.add_argument("input")
    pl.set_defaults(func=cmd_list)

    pa = sub.add_parser("add-primer", help="add a primer_bind annotation")
    pa.add_argument("input")
    pa.add_argument("output")
    pa.add_argument("--name", required=True, help="annotation label")
    pa.add_argument("--seq", required=True, help="primer sequence, 5'->3'")
    pa.add_argument("--strand", choices=["forward", "reverse"], default=None,
                    help="force a strand (default: auto-detect)")
    pa.add_argument("--min-anchor", type=int, default=8,
                    help="minimum perfectly-annealing 3' length (default 8)")
    pa.add_argument("--created-by", default="manual")
    pa.add_argument("--inner-name", default=None,
                    help="inner zip entry name (default: output basename)")
    pa.set_defaults(func=cmd_add_primer)

    pb = sub.add_parser("batch", help="annotate many primers from a CSV/TSV table")
    pb.add_argument("input", help="input .geneious (one sequence)")
    pb.add_argument("output", help="output .geneious (input is never modified)")
    pb.add_argument("primers", help="CSV/TSV with columns: name, sequence")
    pb.add_argument("--min-anchor", type=int, default=15,
                    help="minimum perfect 3' anneal to accept a match (default 15)")
    pb.add_argument("--strand", choices=["forward", "reverse"], default=None,
                    help="restrict to one strand (default: search both)")
    pb.add_argument("--created-by", default="manual")
    pb.add_argument("--inner-name", default=None,
                    help="inner zip entry name (default: output basename)")
    pb.set_defaults(func=cmd_batch)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
