"""Build the final R-vs-S SNV table (candidates + new differential sites + long-read checks).

    python scripts/final_table.py

Inputs: results/summary_candidates.tsv, results/summary_new_diffs_annotated.tsv
(run summarize_batch.py and annotate_diffs.py first). Curated interpretation
from the long-read checks (results/<id>/long/*.tsv) is kept in CURATED below.
Output: results/final_table.tsv
"""
import csv
import os
from pathlib import Path

PROJ = Path(os.environ.get("PROJ", Path(__file__).resolve().parent.parent))
RESULTS = PROJ / "results"
FIXED_AF = 0.8   # same cut-off as af_table.classify_genome_site
PRESENT_AF = 0.05

AMPC = "circuito AmpC / riciclo PG"
# (pair, variant) -> overrides; variant is the aa change or contig:pos for non-coding sites
CURATED = {
    ("254", "S239G"): {
        "af_R": "0.012 (0.005-0.027)", "af_S": "0.991 (0.977-0.996)",
        "esito": "confermato (long reads)", "pathway": "beta-lattamasi",
        "note": "gene in 3 copie; R = 2x CTX-M-15 + 1x CTX-M-178 (copia cromosomica 2.08 Mb, "
                "Ser239), S = 3x CTX-M-15; AF da long reads MAPQ>=20 sulla copia 1"},
    ("246", "N312K"): {"pathway": AMPC,
                       "note": "long reads: esclusiva con P175L (1/418 reads con entrambe); ~62% di R"},
    ("246", "P175L"): {"gene": "dacB (PBP4)", "pathway": AMPC,
                       "note": "seconda mutazione dacB, in sottopopolazione distinta da N312K (~19% di R)"},
    ("216", "G29S"): {"gene": "mpl", "pathway": AMPC,
                      "note": "esclusiva con il frameshift mpl (short reads 63/62/0/0): due sottopopolazioni"},
    ("216", "-1 bp"): {"gene": "mpl", "pathway": AMPC,
                       "note": "in omopolimero G, ma sostenuta dalla fasatura con G29S"},
    ("240", "H190Y"): {"gene": "mpl", "pathway": AMPC},
    ("226", "K136N"): {"pathway": "regolazione porine (EnvZ/OmpR)"},
    ("202", "A97V"): {"pathway": AMPC}, ("240", "H93R"): {"pathway": AMPC},
    ("80", "D135A"): {"pathway": AMPC}, ("236", "H135D"): {"pathway": AMPC},
    ("216", "M1I"): {"pathway": AMPC},
    ("279", "*127E"): {"pathway": "efflusso (MexAB-OprM)"},
    ("226", "L92R"): {"pathway": "involucro / LPS-ECA"},
    ("275", "P952Q"): {"pathway": "RNA polimerasi"},
    ("254", "+1 bp"): {"esito": "probabile artefatto", "note": "omopolimero G, S con 23 reads"},
    ("275", "contig_2:1664570"): {"esito": "probabile artefatto",
                                  "note": "S 57 reads, forte bias di strand (19/2)"},
    ("202", "L466M"): {"esito": "debole", "note": "AF 8%, solo suggestivo"},
    ("240", "V221A"): {"esito": "debole", "note": "AF 8%, solo suggestivo"},
}
OUT_COLS = ["pair", "origine", "gene", "prodotto", "variante", "effetto", "af_R", "af_S",
            "stato", "esito", "pathway", "note"]


def read(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def af_ci(af: str, lo: str, hi: str) -> str:
    return "NA" if af in ("", "NA") else f"{float(af):.3f} ({float(lo):.3f}-{float(hi):.3f})"


def state(af_r: str, af_s: str) -> str:
    vals = {arm: float(v.split()[0]) for arm, v in (("R", af_r), ("S", af_s)) if v != "NA"}
    if len(vals) < 2:
        return "NA"
    arm = max(vals, key=vals.get)
    if vals[arm] < PRESENT_AF:
        return "assente"
    return ("fissata in " if vals[arm] >= FIXED_AF else "subclonale in ") + arm


def candidate_row(c: dict) -> dict:
    af_r, af_s = af_ci(c["af_R"], c["lo_R"], c["hi_R"]), af_ci(c["af_S"], c["lo_S"], c["hi_S"])
    esito = c["verdict"] if c["verdict"] != "non confermato" else "non confermato (artefatto ONT)"
    return {"pair": c["pair"], "origine": "candidato long-read", "gene": c["gene"],
            "prodotto": "", "variante": c["aa"], "effetto": "", "af_R": af_r, "af_S": af_s,
            "esito": esito, "pathway": "", "note": ""}


def new_site_row(d: dict) -> dict:
    variant = d["aa_change"] or f"{d['contig']}:{d['pos']}"
    note = d["context"] if d["feature_type"] == "intergenic" else ""
    return {"pair": d["pair"], "origine": "nuovo (scan short-read)",
            "gene": d["gene"] or d["locus_tag"] or "intergenico", "prodotto": d["product"].replace("%2C", ","),
            "variante": variant, "effetto": d["effect"],
            "af_R": af_ci(d["af_R"], d["lo_R"], d["hi_R"]),
            "af_S": af_ci(d["af_S"], d["lo_S"], d["hi_S"]),
            "esito": "confermato (short reads)", "pathway": "", "note": note}


def finalize(row: dict) -> dict:
    merged = row | CURATED.get((row["pair"], row["variante"]), {})
    return merged | {"stato": state(merged["af_R"], merged["af_S"])}


def main() -> None:
    rows = [candidate_row(c) for c in read(RESULTS / "summary_candidates.tsv")]
    rows += [new_site_row(d) for d in read(RESULTS / "summary_new_diffs_annotated.tsv")]
    final = sorted((finalize(r) for r in rows), key=lambda r: int(r["pair"]))
    out = RESULTS / "final_table.tsv"
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        writer.writeheader()
        writer.writerows(final)
    print(f"[final] {len(final)} rows -> {out}")


if __name__ == "__main__":
    main()
