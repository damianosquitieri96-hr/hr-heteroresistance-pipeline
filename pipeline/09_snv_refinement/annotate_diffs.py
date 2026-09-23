"""Annotate genome-wide differential sites with the bakta annotation of the reference arm.

    python scripts/annotate_diffs.py [results/summary_new_diffs.tsv]

Sites are on the POLISHED reference; bakta annotated the ORIGINAL assembly
(HR_WORK/<id>/annot_<arm>/<arm>.gff3). Each site is lifted back with the same
flank search used for candidates, then placed on the overlapping feature:
CDS -> codon effect (synonymous/missense/nonsense/start_lost/stop_lost,
frameshift/inframe indel); intergenic -> flanking genes, flagged 'promoter'
when within PROMOTER_BP upstream of a gene start.
Output: <input stem>_annotated.tsv next to the input.
"""
import csv
import os
import sys
from pathlib import Path

from liftover_sites import lift, read_fasta

PROJ = Path(os.environ.get("PROJ", Path(__file__).resolve().parent.parent))
SSD = Path(os.environ.get("SSD", PROJ.parent))
PROMOTER_BP = 250
SKIP_TYPES = {"region", "gene", "sequence_feature"}
START_CODONS = {"ATG", "GTG", "TTG"}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: AMINO[16 * i + 4 * j + k]
               for i, a in enumerate(BASES) for j, b in enumerate(BASES)
               for k, c in enumerate(BASES)}
OUT_COLS = ["orig_contig", "orig_pos", "lift_status", "ref_check", "feature_type",
            "locus_tag", "gene", "product", "effect", "aa_change", "context"]


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def parse_attrs(field: str) -> dict:
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


def load_features(gff: Path) -> dict:
    """contig -> list of feature dicts sorted by start (1-based, inclusive)."""
    feats = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] in SKIP_TYPES:
                continue
            a = parse_attrs(f[8])
            feats.setdefault(f[0], []).append({
                "type": f[2], "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                "locus_tag": a.get("locus_tag", a.get("ID", "")),
                "gene": a.get("gene", ""), "product": a.get("product", a.get("Name", ""))})
    return {c: sorted(v, key=lambda x: x["start"]) for c, v in feats.items()}


def cds_effect(seq: str, feat: dict, pos: int, ref: str, alt: str) -> tuple[str, str]:
    """Effect of ref->alt at 1-based pos inside CDS feat (seq = contig sequence)."""
    if len(ref) != len(alt):
        shift = len(alt) - len(ref)
        return ("frameshift" if shift % 3 else "inframe_indel"), f"{shift:+d} bp"
    if len(ref) != 1:
        return "mnp", ""
    cds = seq[feat["start"] - 1:feat["end"]]
    offset = pos - feat["start"]
    if feat["strand"] == "-":
        cds, offset, alt = revcomp(cds), feat["end"] - pos, alt.translate(COMPLEMENT)
    codon_i = offset // 3
    ref_codon = cds[3 * codon_i:3 * codon_i + 3]
    alt_codon = ref_codon[:offset % 3] + alt + ref_codon[offset % 3 + 1:]
    ref_aa, alt_aa = CODON_TABLE.get(ref_codon, "X"), CODON_TABLE.get(alt_codon, "X")
    if codon_i == 0 and ref_codon in START_CODONS:
        ref_aa = "M"
        effect = "synonymous" if alt_codon in START_CODONS else "start_lost"
        alt_aa = "M" if alt_codon in START_CODONS else alt_aa
    elif ref_aa == alt_aa:
        effect = "synonymous"
    elif alt_aa == "*":
        effect = "nonsense"
    elif ref_aa == "*":
        effect = "stop_lost"
    else:
        effect = "missense"
    return effect, f"{ref_aa}{codon_i + 1}{alt_aa}"


def describe_flank(feat: dict, pos: int) -> str:
    gene_start = feat["start"] if feat["strand"] == "+" else feat["end"]
    upstream = (pos < gene_start) if feat["strand"] == "+" else (pos > gene_start)
    dist = abs(pos - gene_start)
    tag = "promoter " if upstream and dist <= PROMOTER_BP else ""
    side = "upstream" if upstream else "downstream"
    name = feat["gene"] or feat["product"]
    return f"{tag}{dist} bp {side} of {feat['locus_tag']} {name} ({feat['strand']})"


def intergenic_context(feats: list, pos: int) -> str:
    left = [f for f in feats if f["end"] < pos]
    right = [f for f in feats if f["start"] > pos]
    parts = [describe_flank(f, pos) for f in (left[-1:] + right[:1])]
    return " | ".join(parts) or "no flanking feature on contig"


def annotate_site(site: dict, orig: dict, polished: dict, feats: dict) -> dict:
    pos, ref, alt = int(site["pos"]), site["ref"].upper(), site["alt"].upper()
    contig, opos, status = lift(polished[site["contig"]], orig, pos)
    row = dict.fromkeys(OUT_COLS, "NA") | {"lift_status": status}
    if contig is None:
        return row
    seq = orig[contig]
    ref_check = "ok" if seq[opos - 1:opos - 1 + len(ref)] == ref else "mismatch"
    row |= {"orig_contig": contig, "orig_pos": opos, "ref_check": ref_check}
    hits = [f for f in feats.get(contig, []) if f["start"] <= opos <= f["end"]]
    if not hits:
        return row | {"feature_type": "intergenic", "effect": "intergenic",
                      "aa_change": "", "locus_tag": "", "gene": "", "product": "",
                      "context": intergenic_context(feats.get(contig, []), opos)}
    feat = next((f for f in hits if f["type"] == "CDS"), hits[0])
    effect, aa = (cds_effect(seq, feat, opos, ref, alt) if feat["type"] == "CDS"
                  else ("non_coding", ""))
    return row | {"feature_type": feat["type"], "locus_tag": feat["locus_tag"],
                  "gene": feat["gene"], "product": feat["product"],
                  "effect": effect, "aa_change": aa, "context": ""}


def pair_inputs(pair: str) -> tuple[dict, dict, dict]:
    ref_dir = PROJ / "results" / pair / "ref"
    origs = list(ref_dir.glob("*.orig.fasta"))
    if len(origs) != 1:
        sys.exit(f"[annotate] expected one *.orig.fasta in {ref_dir}, found {len(origs)}")
    arm = origs[0].name.removesuffix(".orig.fasta")
    gff = SSD / "HR_WORK" / pair / f"annot_{arm}" / f"{arm}.gff3"
    if not gff.is_file():
        sys.exit(f"[annotate] missing annotation {gff}")
    return (read_fasta(origs[0]), read_fasta(ref_dir / f"{arm}.polished.fasta"),
            load_features(gff))


def main() -> None:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else PROJ / "results" / "summary_new_diffs.tsv"
    with open(src) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        in_cols, sites = reader.fieldnames, list(reader)
    cache, rows = {}, []
    for site in sites:
        pair = site["pair"]
        if pair not in cache:
            cache[pair] = pair_inputs(pair)
        rows.append(site | annotate_site(site, *cache[pair]))
    out = src.with_name(src.stem + "_annotated.tsv")
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=in_cols + OUT_COLS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    bad = [r for r in rows if r["ref_check"] != "ok"]
    print(f"[annotate] {len(rows)} sites -> {out}; "
          f"{len(bad)} unlifted or REF mismatch")


if __name__ == "__main__":
    main()
