#!/usr/bin/env python3
"""Transfer candidate SNV coordinates from the original (unpolished) reference
to the polished one.

Polishing can insert/delete bases, so positions shift. For each candidate we
take a flank of the ORIGINAL sequence immediately left of the site and look for
a unique exact match in the polished sequence; if the left flank is not unique
or not found, we try the right flank. The reference base at the lifted position
is reported so the caller can check it still matches the candidate REF.

Output TSV: pair, orig_contig, orig_pos, ref, alt, gene, aa, new_contig,
new_pos, polished_ref_base, lift_status
"""
import argparse
import csv
import sys

FLANK = 60


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    if not seqs:
        sys.exit(f"[liftover] empty FASTA: {path}")
    return seqs


def unique_hits(polished, flank):
    hits = []
    for contig, seq in polished.items():
        start = seq.find(flank)
        while start != -1:
            hits.append((contig, start))
            if len(hits) > 1:
                return hits
            start = seq.find(flank, start + 1)
    return hits


def lift(orig_seq, polished, pos):
    """pos is 1-based on orig_seq. Returns (contig, new_pos_1based, status)."""
    i = pos - 1
    left = orig_seq[max(0, i - FLANK):i]
    if len(left) == FLANK:
        hits = unique_hits(polished, left)
        if len(hits) == 1:
            contig, start = hits[0]
            return contig, start + FLANK + 1, "left_flank"
    right = orig_seq[i + 1:i + 1 + FLANK]
    if len(right) == FLANK:
        hits = unique_hits(polished, right)
        if len(hits) == 1:
            contig, start = hits[0]
            return contig, start, "right_flank"
    return None, None, "unlifted"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidates", required=True, help="semicolon CSV")
    ap.add_argument("--pair", required=True)
    ap.add_argument("--orig", required=True)
    ap.add_argument("--polished", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    orig = read_fasta(args.orig)
    polished = read_fasta(args.polished)

    with open(args.candidates, newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter=";")
                if r["pair"].strip() == args.pair]
    if not rows:
        sys.exit(f"[liftover] no candidates for pair {args.pair}")

    cols = ["pair", "orig_contig", "orig_pos", "ref", "alt", "gene", "aa",
            "prev_afR", "prev_afS", "new_contig", "new_pos",
            "polished_ref_base", "lift_status"]
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            contig, pos = r["contig"].strip(), int(r["pos"])
            if contig not in orig:
                sys.exit(f"[liftover] contig {contig} not in {args.orig}")
            new_contig, new_pos, status = lift(orig[contig], polished, pos)
            base = polished[new_contig][new_pos - 1] if new_contig else "NA"
            if new_contig and base != r["ref"].strip().upper():
                status += ";REF_MISMATCH"
            out.write("\t".join(map(str, [
                args.pair, contig, pos, r["ref"], r["alt"], r["disp"],
                r["aa"], r["afR"], r["afS"], new_contig or "NA",
                new_pos or "NA", base, status])) + "\n")
    print(f"[liftover] {len(rows)} candidates written to {args.out}")


if __name__ == "__main__":
    main()
