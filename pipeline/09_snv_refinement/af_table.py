#!/usr/bin/env python3
"""Allele-frequency tables from a two-sample bcftools mpileup VCF
(FORMAT AD, ADF, ADR), short reads of R and S on the same reference.

No genotype is used: AF = alt reads / (ref + alt reads) per arm, with a Wilson
95% CI, so subclonal variants keep their fraction.

Modes
  candidates  one row per lifted candidate, AF in both arms even when 0,
              plus a verdict against the previous long-read call.
  genome      genome-wide R-vs-S differential sites (new, unbiased scan).
"""
import argparse
import gzip
import math
import sys

from scipy.stats import fisher_exact

MIN_DEPTH = 20          # below this in either arm a site is not assessable
MIN_ALT_READS = 5       # genome scan: alt reads needed in the enriched arm
MIN_ALT_PER_STRAND = 2  # alt must be seen on both strands (DNBSEQ artefact filter)
MIN_DELTA_AF = 0.05     # smallest R-S difference we call differential
P_CANDIDATE = 1e-3      # few pre-specified tests: nominal threshold
P_GENOME = 1e-8         # ~5e6 sites scanned: Bonferroni-like threshold
P_SUGGESTIVE = 1e-5     # genome scan: weaker tier, to be checked by hand
MIN_AF_SUGGESTIVE = 0.03  # enriched-arm AF floor for the suggestive tier
AF_SHIFT_NOTE = 0.2     # |new dAF - old dAF| above this flags a changed fraction
Z95 = 1.96


def wilson(k, n):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + Z95 ** 2 / n
    centre = (p + Z95 ** 2 / (2 * n)) / den
    half = Z95 * math.sqrt(p * (1 - p) / n + Z95 ** 2 / (4 * n * n)) / den
    return (max(0.0, centre - half), min(1.0, centre + half))


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_vcf(path, r_sample, s_sample):
    """Yield (contig, pos, ref, alts, {arm: per-allele dict})."""
    idx = None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                names = f[9:]
                missing = {r_sample, s_sample} - set(names)
                if missing:
                    sys.exit(f"[af_table] samples {missing} not in VCF ({names})")
                idx = {"R": 9 + names.index(r_sample), "S": 9 + names.index(s_sample)}
                continue
            keys = f[8].split(":")
            arms = {}
            for arm, col in idx.items():
                vals = dict(zip(keys, f[col].split(":")))
                arms[arm] = {k: [int(x) if x != "." else 0 for x in vals[k].split(",")]
                             for k in ("AD", "ADF", "ADR")}
            yield f[0], int(f[1]), f[3], f[4].split(","), arms


def arm_stats(a, allele):
    ref, alt = a["AD"][0], a["AD"][allele]
    dp = ref + alt
    lo, hi = wilson(alt, dp)
    return {"dp": dp, "alt": alt, "af": alt / dp if dp else float("nan"),
            "lo": lo, "hi": hi, "altF": a["ADF"][allele], "altR": a["ADR"][allele]}


def compare(arms, allele):
    r, s = arm_stats(arms["R"], allele), arm_stats(arms["S"], allele)
    _, p = fisher_exact([[r["alt"], r["dp"] - r["alt"]], [s["alt"], s["dp"] - s["alt"]]])
    return r, s, p


def fmt(x):
    return f"{x:.4f}" if isinstance(x, float) else str(x)


ARM_COLS = ["dp", "alt", "af", "lo", "hi", "altF", "altR"]


def arm_header(arm):
    return [f"{c}_{arm}" for c in ARM_COLS]


def strand_ok(st):
    return st["altF"] >= MIN_ALT_PER_STRAND and st["altR"] >= MIN_ALT_PER_STRAND


def verdict(r, s, p, prev_r, prev_s):
    if r["dp"] < MIN_DEPTH or s["dp"] < MIN_DEPTH:
        return "non valutabile (profondita')"
    d_new = r["af"] - s["af"]
    enriched = r if d_new > 0 else s
    if p >= P_CANDIDATE or abs(d_new) < MIN_DELTA_AF:
        return "non confermato"
    if not strand_ok(enriched):
        return "non confermato (bias di filamento)"
    d_old = prev_r - prev_s
    if d_old * d_new < 0:
        return "discordante (direzione opposta)"
    if abs(d_new - d_old) > AF_SHIFT_NOTE:
        return "confermato, AF diversa"
    return "confermato"


def load_freebayes(path):
    """(contig, start, end) spans of freebayes records, for overlap lookup.
    freebayes may merge adjacent SNVs into one MNP, so match by overlap."""
    spans = []
    if not path:
        return spans
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            start = int(f[1])
            spans.append((f[0], start, start + len(f[3]) - 1))
    return spans


def fb_hit(spans, contig, pos):
    return any(c == contig and s <= pos <= e for c, s, e in spans)


def run_candidates(args):
    fb = load_freebayes(args.freebayes)
    by_site = {}
    for contig, pos, ref, alts, arms in parse_vcf(args.vcf, args.r_sample, args.s_sample):
        by_site[(contig, pos)] = (ref, alts, arms)

    with open(args.lifted) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cands = [dict(zip(header, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]

    out_cols = header + arm_header("R") + arm_header("S") \
        + ["fisher_p", "dAF_new", "freebayes", "verdict"]
    with open(args.out, "w") as out:
        out.write("\t".join(out_cols) + "\n")
        for c in cands:
            key = (c["new_contig"], int(c["new_pos"])) if c["new_pos"] != "NA" else None
            site = by_site.get(key)
            if site is None:
                out.write("\t".join([c[h] for h in header] + ["NA"] * 17
                                    + ["non valutabile (sito assente)"]) + "\n")
                continue
            fb_called = "si" if fb_hit(fb, *key) else "no"
            ref, alts, arms = site
            alt = c["alt"].strip().upper()
            allele = alts.index(alt) + 1 if alt in alts else None
            if allele is None:
                # alt never observed: AF is 0 in both arms, depth from ref reads
                arms = {k: {"AD": v["AD"][:1] + [0], "ADF": v["ADF"][:1] + [0],
                            "ADR": v["ADR"][:1] + [0]} for k, v in arms.items()}
                allele = 1
            r, s, p = compare(arms, allele)
            v = verdict(r, s, p, float(c["prev_afR"]), float(c["prev_afS"]))
            row = [c[h] for h in header] + [fmt(r[k]) for k in ARM_COLS] \
                + [fmt(s[k]) for k in ARM_COLS] \
                + [f"{p:.3g}", fmt(r["af"] - s["af"]), fb_called, v]
            out.write("\t".join(row) + "\n")
    print(f"[af_table] {len(cands)} candidates -> {args.out}")


def classify_genome_site(r, s, p, d):
    """fixed / subclonal at the genome-wide threshold, subclonal_suggestive at
    the weaker one, None if the site is not differential."""
    enriched = r if d > 0 else s
    if not strand_ok(enriched):
        return None
    if p < P_GENOME and abs(d) >= MIN_DELTA_AF:
        return "fixed" if abs(d) > 0.8 else "subclonal"
    if p < P_SUGGESTIVE and enriched["af"] >= MIN_AF_SUGGESTIVE:
        return "subclonal_suggestive"
    return None


def run_genome(args):
    cols = ["contig", "pos", "ref", "alt"] + arm_header("R") + arm_header("S") \
        + ["fisher_p", "dAF", "class"]
    n_diff, n_fixed_both = {}, 0
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for contig, pos, ref, alts, arms in parse_vcf(args.vcf, args.r_sample, args.s_sample):
            for i, alt in enumerate(alts, start=1):
                if alt == "<*>":
                    continue
                r0, s0 = arm_stats(arms["R"], i), arm_stats(arms["S"], i)
                if max(r0["alt"], s0["alt"]) < MIN_ALT_READS:
                    continue
                if min(r0["dp"], s0["dp"]) < MIN_DEPTH:
                    continue
                if r0["af"] > 0.8 and s0["af"] > 0.8:
                    n_fixed_both += 1   # residual reference error, same in both arms
                    continue
                r, s, p = compare(arms, i)
                d = r["af"] - s["af"]
                cls = classify_genome_site(r, s, p, d)
                if cls is None:
                    continue
                n_diff[cls] = n_diff.get(cls, 0) + 1
                out.write("\t".join([contig, str(pos), ref, alt]
                                    + [fmt(r[k]) for k in ARM_COLS]
                                    + [fmt(s[k]) for k in ARM_COLS]
                                    + [f"{p:.3g}", fmt(d), cls]) + "\n")
    counts = ", ".join(f"{k} {v}" for k, v in sorted(n_diff.items())) or "none"
    print(f"[af_table] genome differential sites: {counts}; "
          f"{n_fixed_both} alt fixed in both arms (residual reference errors)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=["candidates", "genome"])
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--r-sample", required=True)
    ap.add_argument("--s-sample", required=True)
    ap.add_argument("--lifted", help="candidates mode: output of liftover_sites.py")
    ap.add_argument("--freebayes", help="candidates mode: freebayes VCF (cross-check)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if args.mode == "candidates":
        if not args.lifted:
            sys.exit("[af_table] --lifted is required in candidates mode")
        run_candidates(args)
    else:
        run_genome(args)


if __name__ == "__main__":
    main()
