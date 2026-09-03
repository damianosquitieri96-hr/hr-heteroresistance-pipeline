
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory as btf, ScaledTranslation
from matplotlib.patches import Patch
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
try: apply_figure_style(frame="open", sizes=(9,8,7))
except NameError: mpl.rcParams.update({"font.size":9,"savefig.dpi":300})

ip = pd.read_csv(host.artifact_path("3bb29432-83ad-4f3a-9208-e43284c7e73c"))
pw = pd.read_csv(host.artifact_path("2df004f6-0bc3-4d79-a8ca-61591d210020"))
db = pd.read_csv(host.artifact_path("abf66b8a-e419-4777-87ba-24c3d55e435e"))
SP2G = {"E. coli":"Escherichia","K. pneumoniae":"Klebsiella","K. aerogenes":"Klebsiella",
        "E. cloacae":"Enterobacter","P. aeruginosa":"Pseudomonas","C. koseri":"Citrobacter"}
db["genus"] = db.sp2.map(SP2G); assert db.genus.notna().all()
labels = [f"{p}{a}" for p in db.p.astype(str) for a in ("R","S")]
D = pd.DataFrame(0.0, index=labels, columns=labels)
for a,b,d in zip(pw.genome_a, pw.genome_b, pw.mash_dist):
    if a in D.index and b in D.index: D.loc[a,b] = D.loc[b,a] = d
assert (D.values[np.triu_indices(len(labels),1)] > 0).all()
g2 = {f"{p}{a}": g for p,g in zip(db.p.astype(str), db.genus) for a in ("R","S")}
c2 = {f"{p}{a}": ("UKE" if c=="UKE" else "Gemelli") for p,c in zip(db.p.astype(str), db.centre) for a in ("R","S")}
med = float(np.median(ip.loc[ip.pair.astype(str).isin(db.p.astype(str)), "mash"]))
GEN  = {"Escherichia":"#3b6fb0","Klebsiella":"#f4894a","Citrobacter":"#8c7b52",
        "Enterobacter":"#8c66b0","Pseudomonas":"#2aa198"}
SITE = {"UKE":"#5aa457","Gemelli":"#cf5049"}
SPOF = {g: ", ".join(sorted({s for s,gg in SP2G.items() if gg==g})) for g in GEN}
INK, GREY = "#1a1a1a", "#4d4d4d"

ent = [l for l in labels if g2[l] != "Pseudomonas"]
pse = [l for l in labels if g2[l] == "Pseudomonas"]
fig = plt.figure(figsize=(15.6, 12.2))
gs  = fig.add_gridspec(1, 2, width_ratios=[1, 0.60], left=0.235, right=0.955,
                       top=0.915, bottom=0.055, wspace=0.07)
lg  = gs[0].subgridspec(2, 1, height_ratios=[len(ent), len(pse)], hspace=0.115)
axE, axP = fig.add_subplot(lg[0]), fig.add_subplot(lg[1])
rg  = gs[1].subgridspec(4, 2, width_ratios=[1, 0.05], height_ratios=[0.55, 1, 1, 0.55],
                        hspace=0.30, wspace=0.05)
axL = fig.add_subplot(rg[0, :]); axL.axis("off")
axB, axC = fig.add_subplot(rg[1, 0]), fig.add_subplot(rg[2, 0])
cax = fig.add_subplot(rg[1:3, 1]); axN = fig.add_subplot(rg[3, :]); axN.axis("off")

def tree(ax, mem, xmax, ticks, title, letter, xlab):
    sub = D.loc[mem, mem]
    dd = dendrogram(linkage(squareform(sub.values, checks=False), method="average"),
                    labels=mem, no_plot=True)
    ivl = dd["ivl"]; yp = {l: 5+10*i for i, l in enumerate(ivl)}
    for xs, ys in zip(dd["dcoord"], dd["icoord"]): ax.plot(xs, ys, color=GREY, lw=0.8, clip_on=False)
    ax.set_xscale("symlog", linthresh=5e-5, linscale=0.35)
    ax.set_xlim(0, xmax); ax.set_ylim(0, 10*len(ivl))
    ax.set_xticks(ticks)
    ax.set_xticklabels(["0"] + [f"$10^{{{int(np.log10(t))}}}$" for t in ticks[1:]])
    ax.xaxis.set_minor_locator(mpl.ticker.SymmetricalLogLocator(base=10, linthresh=5e-5, subs=range(2,10)))
    ax.set_xlabel(xlab)
    ax.set_yticks([yp[l] for l in ivl]); ax.set_yticklabels(ivl)
    ax.tick_params(axis="y", length=0, pad=26)
    for t in ax.get_yticklabels(): t.set_color(GEN[g2[t.get_text()]]); t.set_fontsize(7)
    ax.spines["left"].set_visible(False)
    offx = lambda inch: btf(ax.transAxes + ScaledTranslation(inch, 0, fig.dpi_scale_trans), ax.transData)
    for s, col in SITE.items():
        ys = [yp[l] for l in ivl if c2[l] == s]
        if ys: ax.plot([0]*len(ys), ys, ls="none", marker="s", ms=6.5, mfc=col, mec="none",
                       transform=offx(-0.245), clip_on=False)
    ax.set_title(title, loc="left", fontsize=9, pad=8)
    ax.text(0, 1.0, letter, transform=btf(ax.transAxes + ScaledTranslation(-1.55, 0, fig.dpi_scale_trans),
            ax.transAxes), fontweight="bold", fontsize=11, va="bottom")
    return yp, offx

ypE, offE = tree(axE, ent, 0.30, [0, 1e-4, 1e-3, 1e-2, 1e-1],
                 "Enterobacterales — 27 R/S pairs, both sites", "a", "Mash distance (21-mer sketch)")
tree(axP, pse, 0.02, [0, 1e-4, 1e-3, 1e-2],
     "$\\it{P.\\ aeruginosa}$ — 6 R/S pairs, Gemelli only", "b", "Mash distance (21-mer sketch)")

def bracket(mem, txt):
    ys = sorted(ypE[m] for m in mem); assert max(ys)-min(ys) == 10*(len(mem)-1), mem
    axE.plot([0,0], [ys[0]-3.5, ys[-1]+3.5], transform=offE(-0.72), color=INK, lw=1.6, clip_on=False)
    axE.text(0, np.mean(ys), txt, transform=offE(-0.82), ha="right", va="center",
             fontsize=8, color=INK, linespacing=1.35, clip_on=False)
bracket(["70R","70S","76R","76S"], "one clone, two patients\nUKE 70 + 76")
bracket(["226R","226S","283R","283S"], "one clone, two patients\nGemelli 226 + 283")
bracket(["87R","87S","220R","220S"], "closest relatives across sites\n$\\it{E.\\ coli}$ 87 (UKE) / 220 (Gemelli)\n"
        f"$2.3\\times10^{{-3}}$ (~{2.252e-3/med:.0f}$\\times$ within-patient)")

fig.suptitle("Every R/S pair is clonal; two patient pairs share a clone, while the closest "
             "cross-site isolates stay ~32$\\times$ more distant", x=0.10, y=0.972, ha="left", fontsize=11)
fig.text(0.10, 0.949, "33 heteroresistance R/S pairs, 66 genomes — UPGMA on the same k-mer sketch as the "
         "gene-dosage figure. Pair 285 excluded (S arm $\\it{M.\\ morganii}$, R arm $\\it{E.\\ coli}$).",
         fontsize=8, color=GREY)
axL.legend(handles=[Patch(facecolor=GEN[g], label=f"$\\it{{{g}}}$ ({SPOF[g]})") for g in
                    ["Escherichia","Klebsiella","Enterobacter","Citrobacter","Pseudomonas"]],
           loc="upper left", bbox_to_anchor=(0, 1.05), frameon=False, handlelength=0.9, handleheight=0.9,
           fontsize=8, title="strain label colour = genus (strain database)", title_fontsize=8, alignment="left")
axL.add_artist(axL.get_legend())
axL.legend(handles=[Patch(facecolor=SITE[s], label=s) for s in ["UKE","Gemelli"]],
           loc="lower left", bbox_to_anchor=(0, -0.10), frameon=False, handlelength=0.9, handleheight=0.9,
           fontsize=8, ncol=2, title="box colour = site", title_fontsize=8, alignment="left")

VMAX = 90.0
for ax, (ttl, mem, letter) in zip([axB, axC],
        [("UKE clone — patients 70 and 76", ["70R","70S","76R","76S"], "c"),
         ("Gemelli clone — patients 226 and 283", ["226R","226S","283R","283S"], "d")]):
    M = D.loc[mem, mem].values*1e6
    im = ax.imshow(np.ma.masked_array(M, mask=np.eye(4, dtype=bool)), cmap="Oranges", vmin=0, vmax=VMAX)
    im.cmap.set_bad("#ededed")
    for i in range(4):
        for j in range(4):
            ax.text(j, i, "\u2014" if i==j else f"{M[i,j]:.0f}", ha="center", va="center", fontsize=8,
                    color="#8a8a8a" if i==j else ("white" if M[i,j] > 0.62*VMAX else "black"))
    ax.set_xticks(range(4)); ax.set_yticks(range(4))
    ax.set_xticklabels(mem); ax.set_yticklabels(mem); ax.tick_params(length=0)
    ax.set_xticks(np.arange(-.5,4), minor=True); ax.set_yticks(np.arange(-.5,4), minor=True)
    ax.grid(which="minor", color="white", lw=2); ax.tick_params(which="minor", length=0)
    for s in ax.spines.values(): s.set_visible(False)
    for t, c in zip(ax.get_xticklabels()+ax.get_yticklabels(), [SITE[c2[m]] for m in mem]*2): t.set_color(c)
    ax.set_title(ttl, loc="left", fontsize=9, pad=6)
    ax.text(-0.17, 1.13, letter, transform=ax.transAxes, fontweight="bold", fontsize=11)
cb = fig.colorbar(im, cax=cax); cb.outline.set_visible(False)
cb.set_label("Mash distance, $\\times10^{-6}$ (same scale in c and d)", fontsize=8)
cb.set_ticks([0, 30, 60, 90])
cax.axhline(med*1e6, color="black", lw=1.1, ls=(0,(3,2)))

axN.text(0, 1.0, "Gemelli clone: 283R and 283S are each closer to 226S (62, 61)\n"
         "than to one another (66) — the only two genomes whose nearest\n"
         "neighbour is not its own R/S partner (64/66 pass); wells 41/42 and\n"
         "63/64 are not adjacent, so plate carry-over does not explain it.\n"
         "UKE clone: cross-patient distances (80–88) exceed the within-\n"
         "patient ones (57, 48) — each patient carries its own variant.\n"
         "Dashed line on the colour scale = cohort median within-patient\n"
         f"distance ({med*1e6:.0f}).",
         transform=axN.transAxes, ha="left", va="top", fontsize=8, color="#333333", linespacing=1.5)
fig.savefig("HR_clonality_onetree.png", dpi=300, facecolor="white")
r = fig.canvas.get_renderer()
tx = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print("overlaps:", len([1 for i,(a,ba) in enumerate(tx) for b,bb in tx[i+1:] if ba.overlaps(bb)]),
      "| clipped:", [t.get_text()[:22] for t,b in tx if b.x0 < 0 or b.x1 > fig.bbox.x1])
