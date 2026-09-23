"""Long-read checks on the polished reference (BAMs from map_long.sh).

  phase   haplotype counts of 2+ nearby sites on reads spanning all of them
          python scripts/long_read_check.py phase --bam B --site contig:pos:REF:ALT ...
  copies  base counts at one site in each copy of a repeated gene, per arm
          python scripts/long_read_check.py copies --bam B1 --bam B2 --site contig:pos:REF:ALT ...

Primary alignments only, MAPQ >= MAPQ_MIN. A read's call at a site is 'ref',
'alt' or 'other'; indel sites are called from the net CIGAR indel length over
the whole homopolymer/ref span, SNV sites from the aligned base.
"""
import argparse
import re
import subprocess
from collections import Counter
from itertools import product

from af_table import wilson

MAPQ_MIN = 20
SAMTOOLS = "samtools"
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = set("MDN=X")
QUERY_CONSUMING = set("MIS=X")


def parse_site(text: str) -> dict:
    contig, pos, ref, alt = text.rsplit(":", 3)
    return {"contig": contig, "pos": int(pos), "ref": ref.upper(), "alt": alt.upper(),
            "label": text}


def aligned_reads(bam: str, contig: str, start: int, end: int):
    """Yield (name, ref_start_1based, cigar_ops, seq) of primary reads covering start..end."""
    cmd = [SAMTOOLS, "view", "-F", "0x904", "-q", str(MAPQ_MIN), bam, f"{contig}:{start}-{end}"]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        ops = [(int(n), op) for n, op in CIGAR_RE.findall(f[5])]
        ref_len = sum(n for n, op in ops if op in REF_CONSUMING)
        read_start = int(f[3])
        if read_start <= start and read_start + ref_len - 1 >= end:
            yield f[0], read_start, ops, f[9]


def walk(read_start: int, ops: list):
    """Yield (op, length, ref_pos_1based, query_index) for each CIGAR block."""
    rpos, qpos = read_start, 0
    for n, op in ops:
        yield op, n, rpos, qpos
        rpos += n if op in REF_CONSUMING else 0
        qpos += n if op in QUERY_CONSUMING else 0


def base_at(read_start: int, ops: list, seq: str, pos: int) -> str:
    for op, n, rpos, qpos in walk(read_start, ops):
        if op in REF_CONSUMING and rpos <= pos < rpos + n:
            return seq[qpos + pos - rpos] if op in "M=X" else "-"
    return "?"


def net_indel(read_start: int, ops: list, start: int, end: int) -> int:
    """Inserted minus deleted bases with the event anchored within [start, end]."""
    net = 0
    for op, n, rpos, _ in walk(read_start, ops):
        if op == "I" and start <= rpos <= end + 1:
            net += n
        elif op == "D" and rpos <= end and rpos + n - 1 >= start:
            net -= n
    return net


def call(site: dict, read_start: int, ops: list, seq: str) -> str:
    ref, alt = site["ref"], site["alt"]
    if len(ref) == len(alt) == 1:
        b = base_at(read_start, ops, seq, site["pos"])
        return "ref" if b == ref else "alt" if b == alt else "other"
    net = net_indel(read_start, ops, site["pos"], site["pos"] + len(ref) - 1)
    return "ref" if net == 0 else "alt" if net == len(alt) - len(ref) else "other"


def fmt_af(k: int, n: int) -> str:
    lo, hi = wilson(k, n)
    return f"{k}/{n} = {k / n:.3f} ({lo:.3f}-{hi:.3f})" if n else "0/0"


def run_phase(args) -> None:
    sites = [parse_site(s) for s in args.site]
    contig = {s["contig"] for s in sites}
    if len(contig) != 1:
        raise SystemExit("[phase] all sites must be on one contig")
    start = min(s["pos"] for s in sites)
    end = max(s["pos"] + len(s["ref"]) - 1 for s in sites)
    haps = Counter()
    for _, rs, ops, seq in aligned_reads(args.bam, contig.pop(), start, end):
        haps[tuple(call(s, rs, ops, seq) for s in sites)] += 1
    total = sum(haps.values())
    print(f"# {args.bam}\n# reads spanning all sites (MAPQ>={MAPQ_MIN}): {total}")
    for i, s in enumerate(sites):
        k = sum(v for h, v in haps.items() if h[i] == "alt")
        n = sum(v for h, v in haps.items() if h[i] in ("ref", "alt"))
        print(f"site {s['label']}: long-read AF {fmt_af(k, n)}; other={total - n}")
    print("haplotype\t" + "\t".join(s["label"] for s in sites) + "\treads")
    for hap in product(("ref", "alt", "other"), repeat=len(sites)):
        if haps[hap]:
            print("\t" + "\t".join(hap) + f"\t{haps[hap]}")


def run_copies(args) -> None:
    sites = [parse_site(s) for s in args.site]
    print("bam\tcopy\treads_mapq\tref\talt\tother\talt_AF")
    for bam in args.bam:
        for s in sites:
            counts = Counter(call(s, rs, ops, seq)
                             for _, rs, ops, seq in aligned_reads(bam, s["contig"], s["pos"], s["pos"]))
            n = counts["ref"] + counts["alt"]
            print(f"{bam.rsplit('/', 1)[-1]}\t{s['label']}\t{sum(counts.values())}\t"
                  f"{counts['ref']}\t{counts['alt']}\t{counts['other']}\t{fmt_af(counts['alt'], n)}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("phase")
    p.add_argument("--bam", required=True)
    p.add_argument("--site", action="append", required=True)
    c = sub.add_parser("copies")
    c.add_argument("--bam", action="append", required=True)
    c.add_argument("--site", action="append", required=True)
    args = ap.parse_args()
    run_phase(args) if args.cmd == "phase" else run_copies(args)


if __name__ == "__main__":
    main()
