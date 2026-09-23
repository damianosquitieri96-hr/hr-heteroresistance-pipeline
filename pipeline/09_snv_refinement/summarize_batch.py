"""Merge per-pair refinement outputs into one summary table.

    python scripts/summarize_batch.py

Reads results/<id>/candidates_AF.tsv and results/<id>/genome_differential.tsv
for every pair folder and writes:
  results/summary_candidates.tsv   one row per original candidate (AF R/S, verdict)
  results/summary_new_diffs.tsv    genome-wide differential sites NOT among the candidates
  results/summary_pairs.tsv        per-pair counts
"""
import os
import sys
from pathlib import Path

import pandas as pd

PROJ = Path(os.environ.get("PROJ", Path(__file__).resolve().parent.parent))
RESULTS = PROJ / "results"
CAND_COLS = ["pair", "gene", "aa", "orig_pos", "ref", "alt", "prev_afR", "prev_afS",
             "dp_R", "af_R", "lo_R", "hi_R", "dp_S", "af_S", "lo_S", "hi_S",
             "fisher_p", "freebayes", "verdict"]


def pair_dirs() -> list[Path]:
    return sorted((d for d in RESULTS.iterdir() if d.is_dir() and d.name.isdigit()),
                  key=lambda d: int(d.name))


def read_tsv(path: Path) -> pd.DataFrame | None:
    if not path.is_file():
        print(f"[warn] missing {path.relative_to(PROJ)}", file=sys.stderr)
        return None
    return pd.read_csv(path, sep="\t")


def new_diffs(pair: str, diff: pd.DataFrame, cand: pd.DataFrame | None) -> pd.DataFrame:
    known = set() if cand is None else set(zip(cand["new_contig"], cand["new_pos"]))
    mask = [(c, p) not in known for c, p in zip(diff["contig"], diff["pos"])]
    return diff.loc[mask].assign(pair=pair)


def main() -> None:
    cand_frames, diff_frames, pair_rows = [], [], []
    for d in pair_dirs():
        cand = read_tsv(d / "candidates_AF.tsv")
        diff = read_tsv(d / "genome_differential.tsv")
        if cand is not None:
            cand_frames.append(cand[CAND_COLS])
        extra = new_diffs(d.name, diff, cand) if diff is not None else pd.DataFrame()
        if not extra.empty:
            diff_frames.append(extra)
        verdicts = cand["verdict"].value_counts().to_dict() if cand is not None else {}
        pair_rows.append({"pair": d.name,
                          "n_candidates": 0 if cand is None else len(cand),
                          **{f"n_{k.replace(' ', '_')}": v for k, v in verdicts.items()},
                          "n_genome_diff": 0 if diff is None else len(diff),
                          "n_new_diff": len(extra),
                          "complete": cand is not None and diff is not None})

    if cand_frames:
        pd.concat(cand_frames).to_csv(RESULTS / "summary_candidates.tsv", sep="\t", index=False)
    if diff_frames:
        cols = ["pair"] + [c for c in diff_frames[0].columns if c != "pair"]
        pd.concat(diff_frames)[cols].to_csv(RESULTS / "summary_new_diffs.tsv", sep="\t", index=False)
    pairs = pd.DataFrame(pair_rows).fillna(0)
    pairs.to_csv(RESULTS / "summary_pairs.tsv", sep="\t", index=False)
    print(pairs.to_string(index=False))


if __name__ == "__main__":
    main()
