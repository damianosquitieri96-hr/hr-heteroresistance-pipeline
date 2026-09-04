"""Defensible intra-pair chromosomal substitution distance, from the joint R-vs-S VCFs.

Rationale (spec163 methodological re-check): counts obtained by assembly-vs-assembly
BLASTN over sliding windows inflate the distance in a non-constant way, because
paralogue misalignment inside repeat families contributes mismatches that move when
the window boundaries move. The exports already contain read-based joint calls
(bcftools, samples R and S) against ONE ARM's own assembly from the same pair, which is
the within-pair reference choice recommended in Gorrie et al. Lancet Microbe
2021;2:e575-83, so the distance is recomputed from those instead of re-running an
assembly comparison.

Rules applied, one per point of the proposal:
 1. reference = an assembly of the same pair. It is the S arm in 8 of the 17 pairs with
    a joint call set and the R arm in the other 9 (data/pair_reference_arm.csv, read
    from the VCF ##reference headers). Nothing below depends on which: samples are
    located by their R/S suffix in the VCF header, and classify() is symmetric in the
    two alt fractions, so a difference is detected with either arm as reference.
    CAVEAT for readers of the output: af_R and af_S are alt fractions relative to THAT
    PAIR'S reference arm, so their polarity flips between pairs. At a fixed difference,
    an S-reference pair shows af_R ~ 0.99 / af_S ~ 0.00 and an R-reference pair shows
    af_R ~ 0.00 / af_S ~ 1.00. Filtering the differential-sites table on af_R >= 0.9 to
    get "variants carried by the resistant arm" would therefore silently drop every
    R-reference pair; use cls together with the reference arm instead.
 2. no prophage and no recombination masking; declared, not applied.
 3. chromosome only: the largest reference contig (any further contig >= 1 Mb is
    reported separately, never silently pooled). Plasmid contigs are excluded.
 4. dense divergent tracts (>= 5 differential sites within 500 bp) are collapsed to
    one event and reported apart from dispersed substitutions.
 5. everything discarded is counted, so the audit trail is complete: shared alleles,
    low-depth sites, indels.

Difference classes, from per-sample allelic depths (FORMAT/AD), both arms >= MINDP:
   fixed      one arm >= 0.90 alt fraction, the other <= 0.10
   subclonal  one arm >= 0.20, the other <= 0.05, and the two differ by >= 0.20
   shared     neither of the above -> not a difference between the arms
"""
import gzip, glob, os
import pandas as pd

E = "/Users/damsq/Downloads/hr_exports"
MINDP, HI, LO, SUB_HI, SUB_LO, DSUB = 10, 0.90, 0.10, 0.20, 0.05, 0.20
TRACT_N, TRACT_BP, CHR_MIN = 5, 500, 1_000_000


def af(sample_field, fmt):
    """alt fraction and depth from a bcftools FORMAT field."""
    keys = fmt.split(":")
    if "AD" not in keys:
        return None, 0
    ad = [int(x) for x in sample_field.split(":")[keys.index("AD")].split(",") if x not in (".", "")]
    dp = sum(ad)
    return (None, 0) if dp == 0 else (1 - ad[0] / dp, dp)


def classify(afr, afs):
    hi, lo = max(afr, afs), min(afr, afs)
    if hi >= HI and lo <= LO:
        return "fixed"
    if hi >= SUB_HI and lo <= SUB_LO and (hi - lo) >= DSUB:
        return "subclonal"
    return "shared"


def recount(pair):
    vcf = glob.glob(f"{E}/{pair}/*.RS.vcf.gz")[0]
    ctg = pd.read_csv(glob.glob(f"{E}/{pair}/*.contigs.tsv")[0], sep="\t", header=None,
                      names=["contig", "length"])
    chrom = ctg.sort_values("length", ascending=False).iloc[0]
    extra = ctg[(ctg.length >= CHR_MIN) & (ctg.contig != chrom.contig)]

    with gzip.open(vcf, "rt") as fh:
        rows, stats = [], dict(low_depth=0, shared=0, indel_diff=0, off_chrom=0, no_ad=0)
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = [os.path.basename(s).split(".")[0] for s in f[9:]]
                iR = next(i for i, s in enumerate(samples) if s.endswith("R"))
                iS = next(i for i, s in enumerate(samples) if s.endswith("S"))
                continue
            if f[0] != chrom.contig:
                stats["off_chrom"] += 1
                continue
            afr, dpr = af(f[9 + iR], f[8])
            afs, dps = af(f[9 + iS], f[8])
            if afr is None or afs is None:
                stats["no_ad"] += 1
                continue
            if dpr < MINDP or dps < MINDP:
                stats["low_depth"] += 1
                continue
            cls = classify(afr, afs)
            if cls == "shared":
                stats["shared"] += 1
                continue
            if not (len(f[3]) == 1 and all(len(a) == 1 for a in f[4].split(","))):
                stats["indel_diff"] += 1
                continue
            rows.append(dict(pair=pair, contig=f[0], pos=int(f[1]), ref=f[3], alt=f[4],
                             af_R=round(afr, 3), af_S=round(afs, 3), dp_R=dpr, dp_S=dps, cls=cls))

    D = pd.DataFrame(rows)
    tracts = 0
    if len(D):
        D = D.sort_values("pos").reset_index(drop=True)
        pos = D.pos.to_numpy()
        grp, g = [0], 0
        for i in range(1, len(pos)):
            if pos[i] - pos[i - 1] > TRACT_BP:
                g += 1
            grp.append(g)
        D["grp"] = grp
        sz = D.groupby("grp").size()
        dense = set(sz[sz >= TRACT_N].index)
        tracts = len(dense)
        D["in_tract"] = D.grp.isin(dense)
    else:
        D = pd.DataFrame(columns=["pair", "contig", "pos", "ref", "alt", "af_R", "af_S",
                                  "dp_R", "dp_S", "cls", "grp", "in_tract"])

    disp = D[~D.in_tract.astype(bool)]
    return D, dict(pair=pair, chrom=chrom.contig, chrom_bp=int(chrom.length),
                   extra_large_contigs=";".join(extra.contig) or "-",
                   fixed_dispersed=int((disp.cls == "fixed").sum()),
                   subclonal_dispersed=int((disp.cls == "subclonal").sum()),
                   dense_tracts=tracts, sites_in_tracts=int(D.in_tract.astype(bool).sum()),
                   indel_differences=stats["indel_diff"], shared_sites=stats["shared"],
                   low_depth_sites=stats["low_depth"], off_chromosome=stats["off_chrom"])


pairs = sorted({os.path.basename(os.path.dirname(v)) for v in glob.glob(f"{E}/*/*.RS.vcf.gz")},
               key=lambda x: int(x))
alldiff, summ = [], []
for p in pairs:
    D, s = recount(p)
    alldiff.append(D)
    s["defensible_distance"] = s["fixed_dispersed"] + s["dense_tracts"]
    summ.append(s)

S = pd.DataFrame(summ)
S.to_csv("snpdist/HR_intrapair_chromosomal_distance.csv", index=False)
pd.concat(alldiff, ignore_index=True).to_csv("snpdist/HR_intrapair_differential_sites.csv", index=False)
print(S.to_string(index=False))
