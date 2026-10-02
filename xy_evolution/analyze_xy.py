#!/usr/bin/env python3
"""Analyse mouse chrX / chrY homology and conservation across vertebrates.

Inputs  (from fetch_xy_data.py / phastcons_bins.py)
  data/cache/mouse_genes.json, data/cache/orthologs.json, data/cache/cds.fna.gz
  data/phastcons_zoom.json, ../sequences/*.genes.tsv (Yq multicopy families)
Outputs
  data/genes.tsv        one row per mouse X/Y protein-coding gene
  data/orthologs.tsv    one row per mouse gene x panel species (best ortholog)
  data/gametologs.tsv   Y gene -> most similar X gene, protein identity, dS, dN
  data/windows.tsv      1 Mb windows: gene density, function, conservation
  data/xy_data.json     everything the viewer needs

Per mouse gene
  * protein identity to the ortholog in each species (global alignment of the
    longest-CDS translations, BLOSUM62, end gaps free; identity over aligned
    residue pairs, coverage = aligned residues / mouse length)
  * dN, dS (Nei-Gojobori 1986 + Jukes-Cantor) to rat and human
  * phylostratum = most distant clade with an ortholog
  * multicopy family = connected component of mouse X/Y proteins >= 70 %
    identical over >= 50 % of the shorter one
  * GO: experimental evidence (EXP/IDA/IPI/IMP/IGI/IEP/HDA/HMP/HGI/HEP) and
    functional classes (reproduction, nervous system, immune, chromatin/RNA)
  * phastCons35way mean over the gene body
"""
import collections
import csv
import glob
import gzip
import json
import math
import multiprocessing as mp
import os
import sys

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch_xy_data import ENS_SPECIES, PANEL, XY_PAIRS  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
CACHE = os.path.join(DATA, "cache")
WINDOW = 1_000_000
CHROM_LEN = {"X": 169_476_592, "Y": 91_455_967}  # GRCm39
STRATA = [  # clade order for phylostratum, closest first
    ("Murinae", ["rat"]), ("Rodentia", ["hamster", "guinea_pig"]),
    ("Euarchontoglires", ["human", "chimp", "macaque", "marmoset"]),
    ("Eutheria", ["dog", "cow", "horse", "elephant"]), ("Theria", ["opossum"]),
    ("Mammalia", ["platypus"]), ("Amniota", ["chicken", "anole"]),
    ("Tetrapoda", ["xenopus"]), ("Vertebrata", ["zebrafish"]),
]
EXPERIMENTAL = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HDA", "HMP", "HGI", "HEP", "HTP"}
FUNC_CLASSES = {
    "reproduction": ["spermat", "sperm", "reproduct", "gamete", "meio", "fertiliz", "male gonad",
                     "testis", "oocyte", "ovarian", "sex determin", "acrosom", "spermatid",
                     "seminiferous", "germ cell", "piRNA", "female pregnancy", "placenta"],
    "nervous": ["neuro", "synap", "brain", "axon", "dendrit", "learning", "memory", "cognition",
                "behavior", "nervous system", "glial", "myelin", "cerebral", "hippocamp"],
    "immune": ["immun", "lymphocyte", "T cell", "B cell", "leukocyte", "inflamm", "cytokine",
               "interferon", "antigen", "defense response"],
    "chromatin_rna": ["chromatin", "histone", "RNA splicing", "mRNA processing", "translation",
                      "transcription by RNA polymerase", "DNA repair", "ubiquitin"],
}

_AL = {}


def aligner():
    if "a" not in _AL:
        a = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                            open_gap_score=-11, extend_gap_score=-1)
        a.end_gap_score = 0
        _AL["a"] = a
    return _AL["a"]


def translate(cds):
    p = str(Seq(cds).translate())
    return p[:-1] if p.endswith("*") else p


# ---------------------------------------------------------------- dN / dS (NG86)
BASES = "TCAG"
CODE = dict(zip((a + b + c for a in BASES for b in BASES for c in BASES),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def syn_sites(c):
    if CODE.get(c, "*") == "*":
        return 0.0
    s = 0.0
    for i in range(3):
        for b in "TCAG":
            if b != c[i]:
                aa = CODE[c[:i] + b + c[i + 1:]]
                if aa == CODE[c]:
                    s += 1 / 3
    return s


SYN_SITES = {c: syn_sites(c) for c in CODE}
_PATH_CACHE = {}


def codon_diffs(c1, c2):
    """(synonymous, nonsynonymous) differences averaged over mutational paths."""
    key = (c1, c2)
    if key in _PATH_CACHE:
        return _PATH_CACHE[key]
    pos = [i for i in range(3) if c1[i] != c2[i]]
    import itertools
    sd = nd = 0.0
    n = 0
    for order in itertools.permutations(pos):
        cur, s, d, ok = c1, 0, 0, True
        for i in order:
            nxt = cur[:i] + c2[i] + cur[i + 1:]
            if CODE[nxt] == "*":
                ok = False
                break
            if CODE[nxt] == CODE[cur]:
                s += 1
            else:
                d += 1
            cur = nxt
        if ok:
            sd, nd, n = sd + s, nd + d, n + 1
    res = (sd / n, nd / n) if n else (0.0, float(len(pos)))
    _PATH_CACHE[key] = res
    return res


def ng86(cds1, cds2, aln):
    """dN, dS from a protein alignment of the two translated CDSs."""
    S = N = Sd = Nd = 0.0
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for k in range(e1 - s1):
            c1 = cds1[3 * (s1 + k):3 * (s1 + k) + 3]
            c2 = cds2[3 * (s2 + k):3 * (s2 + k) + 3]
            if len(c1) < 3 or len(c2) < 3 or c1 not in CODE or c2 not in CODE:
                continue
            if CODE[c1] == "*" or CODE[c2] == "*":
                continue
            ss = (SYN_SITES[c1] + SYN_SITES[c2]) / 2
            S += ss
            N += 3 - ss
            if c1 != c2:
                a, b = codon_diffs(c1, c2)
                Sd += a
                Nd += b
    if S < 10 or N < 10:
        return None, None

    def jc(p):
        x = 1 - 4 * p / 3
        return -0.75 * math.log(x) if x > 0 else None

    return jc(Nd / N), jc(Sd / S)


def compare(args):
    """Align two CDSs (via proteins). Returns identity, coverage, dN, dS."""
    cds1, cds2, want_dnds = args
    p1, p2 = translate(cds1), translate(cds2)
    if not p1 or not p2 or len(p1) * len(p2) > 2.5e7:
        return None
    try:
        aln = aligner().align(p1, p2)[0]
    except OverflowError:
        return None
    match = pairs = 0
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        pairs += e1 - s1
        match += sum(1 for a, b in zip(p1[s1:e1], p2[s2:e2]) if a == b)
    if not pairs:
        return None
    out = {"pid": round(100.0 * int(match) / int(pairs), 1), "cov": round(float(pairs) / len(p1), 3),
           "cov2": round(float(pairs) / len(p2), 3)}
    if want_dnds:
        dn, ds = ng86(cds1, cds2, aln)
        out["dn"] = None if dn is None else abs(round(dn, 4))
        out["ds"] = None if ds is None else abs(round(ds, 4))
    return out


# ---------------------------------------------------------------- helpers
def load_cds():
    cds = {}
    with gzip.open(os.path.join(CACHE, "cds.fna.gz"), "rt") as f:
        gid = None
        for line in f:
            if line.startswith(">"):
                gid = line[1:].strip()
            else:
                cds[gid] = line.strip()
    return cds


def kmers(p, k=5):
    return {p[i:i + k] for i in range(len(p) - k + 1)}


def similar_pairs(prot, min_shared=0.2):
    """Candidate protein pairs sharing >= min_shared of the shorter one's 5-mers."""
    index = collections.defaultdict(set)
    km = {g: kmers(p) for g, p in prot.items()}
    for g, ks in km.items():
        for k in ks:
            index[k].add(g)
    pairs = set()
    for g, ks in km.items():
        cnt = collections.Counter()
        for k in ks:
            if len(index[k]) < 400:
                cnt.update(index[k])
        for h, n in cnt.items():
            if h > g and n >= min_shared * min(len(ks), len(km[h])) and n >= 5:
                pairs.add((g, h))
    return pairs


def func_classes(go):
    names = [t[2] for t in go if t[0] == "B"]
    out = []
    for cls, keys in FUNC_CLASSES.items():
        if any(k.lower() in n.lower() for n in names for k in keys):
            out.append(cls)
    return out


def phast_index(zoom):
    """Prefix sums over zoom records for fast range means."""
    idx = {}
    for chrom in ("chrX", "chrY"):
        recs = zoom[chrom]
        starts = [r[0] for r in recs]
        cn = [0]
        cs = [0.0]
        for r in recs:
            cn.append(cn[-1] + r[2])
            cs.append(cs[-1] + r[3])
        idx[chrom[3:]] = (starts, recs, cn, cs)
    return idx


def phast_mean(idx, chrom, start, end):
    import bisect
    starts, recs, cn, cs = idx[chrom]
    i = max(bisect.bisect_right(starts, start) - 1, 0)
    j = bisect.bisect_left(starts, end)
    n = cn[j] - cn[i]
    return (cs[j] - cs[i]) / n if n > 0 else None


def yq_families():
    fam = {}
    for path in glob.glob(os.path.join(HERE, "..", "sequences", "*.genes.tsv")):
        f = os.path.basename(path).split(".")[0]
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                fam[row["gene_id"]] = f
    return fam


# ---------------------------------------------------------------- main
def main():
    mouse = json.load(open(os.path.join(CACHE, "mouse_genes.json")))
    orth = json.load(open(os.path.join(CACHE, "orthologs.json")))
    cds = load_cds()
    zoom = json.load(open(os.path.join(DATA, "phastcons_zoom.json")))
    pidx = phast_index(zoom)
    yq = yq_families()
    sp_name = {t: v[0] for t, v in PANEL.items()}
    sp_clade = {v[0]: v[1] for v in PANEL.values()}

    # ---- 1. pairwise ortholog comparisons
    jobs, keys = [], []
    for gid, olist in orth.items():
        if gid not in cds:
            continue
        for o in olist:
            if o["gene_id"] in cds:
                sp = sp_name[o["tax_id"]]
                jobs.append((cds[gid], cds[o["gene_id"]], sp in ("rat", "human")))
                keys.append((gid, sp, o))
    print(f"{len(jobs)} ortholog alignments", file=sys.stderr)
    with mp.Pool() as pool:
        res = pool.map(compare, jobs, chunksize=20)
    best = {}
    for (gid, sp, o), r in zip(keys, res):
        rec = {"sym": o["symbol"], "gid": o["gene_id"],
               "chr": (o["loc"] or {}).get("chr") or ",".join(o["chromosomes"]) or "?",
               "start": (o["loc"] or {}).get("start"), "n": 1, **(r or {})}
        prev = best.get((gid, sp))
        if prev is None or rec.get("pid", -1) > prev.get("pid", -1):
            if prev:
                rec["n"] = prev["n"] + 1
            best[(gid, sp)] = rec
        else:
            prev["n"] += 1
    # Ensembl Compara supplement: species where NCBI Orthologs has no ortholog
    ens = json.load(open(os.path.join(CACHE, "ensembl_orthologs.json")))
    ens_loc = json.load(open(os.path.join(CACHE, "ensembl_locations.json")))
    ens_sp = {v: sp_name[t] for t, v in ENS_SPECIES.items()}
    n_added = 0
    for gid, g in mouse.items():
        cands = collections.defaultdict(list)
        for e in g.get("ens", []):
            for h in ens.get(e) or []:
                if h["species"] in ens_sp:
                    cands[ens_sp[h["species"]]].append(h)
        for sp, hs in cands.items():
            if (gid, sp) in best:
                continue
            h = max(hs, key=lambda h: max(h["pid_source"], h["pid_target"]))
            if max(h["pid_source"], h["pid_target"]) < 25:  # distant one2many tree calls
                continue
            loc = ens_loc.get(h["id"], {})
            best[(gid, sp)] = {"sym": loc.get("sym", h["id"]), "gid": h["id"], "chr": loc.get("chr") or "?",
                               "start": loc.get("start"), "n": len(hs), "src": "Ensembl",
                               "pid": round(max(h["pid_source"], h["pid_target"]), 1)}
            n_added += 1
    print(f"{n_added} orthologs added from Ensembl Compara", file=sys.stderr)
    for gid, olist in orth.items():  # orthologs without CDS still count for presence
        for o in olist:
            sp = sp_name[o["tax_id"]]
            if (gid, sp) not in best:
                best[(gid, sp)] = {"sym": o["symbol"], "gid": o["gene_id"],
                                   "chr": (o["loc"] or {}).get("chr") or ",".join(o["chromosomes"]) or "?",
                                   "start": (o["loc"] or {}).get("start"), "n": 1}

    # ---- 2. multicopy families among mouse X/Y proteins
    prot = {g: translate(cds[g]) for g in mouse if g in cds}
    prot = {g: p for g, p in prot.items() if len(p) >= 30}
    cand = sorted(similar_pairs(prot))
    print(f"{len(cand)} paralog candidate pairs", file=sys.stderr)
    with mp.Pool() as pool:
        pres = pool.map(compare, [(cds[a], cds[b], False) for a, b in cand], chunksize=20)
    parent = {g: g for g in prot}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    para_best = collections.defaultdict(lambda: (0, None))
    for (a, b), r in zip(cand, pres):
        if not r:
            continue
        same_chr = (mouse[a]["loc"] or {}).get("chr") == (mouse[b]["loc"] or {}).get("chr")
        if same_chr and r["pid"] >= 70 and max(r["cov"], r["cov2"]) >= 0.5:
            parent[find(a)] = find(b)
        for x, y in ((a, b), (b, a)):
            if r["pid"] * min(r["cov"], r["cov2"]) > para_best[x][0]:
                para_best[x] = (r["pid"] * min(r["cov"], r["cov2"]), y, r)
    comp = collections.defaultdict(list)
    for g in prot:
        comp[find(g)].append(g)

    def fam_name(members):
        syms = [mouse[m]["symbol"] for m in members]
        named = [s for s in syms if not s.startswith(("Gm", "LOC"))]
        for m in members:
            if m in yq:
                return yq[m]
        return sorted(named or syms, key=lambda s: (len(s), s))[0] + "_fam"

    family = {}
    for members in comp.values():
        if len(members) >= 2:
            name = fam_name(members)
            for m in members:
                family[m] = (name, len(members))
    for g, f in yq.items():  # keep the curated Yq assignment (copies without CDS too)
        if g in mouse:
            family[g] = (f, sum(1 for v in yq.values() if v == f))

    # ---- 3. X-Y gametologs: each Y gene vs its most similar X gene
    xprot = {g: p for g, p in prot.items() if (mouse[g]["loc"] or {}).get("chr") == "X"}
    yprot = {g: p for g, p in prot.items() if g not in xprot}
    # looser prefilter than for paralogs: X-Y pairs can be < 40 % identical
    xy_cand = collections.defaultdict(set)
    index = collections.defaultdict(set)
    for g, p in xprot.items():
        for k in kmers(p, 3):
            index[k].add(g)
    xk = {g: len(kmers(p, 3)) for g, p in xprot.items()}
    for y, p in yprot.items():
        cnt = collections.Counter()
        for k in kmers(p, 3):
            if len(index[k]) < 200:
                cnt.update(index[k])
        hits = [(x, n) for x, n in cnt.items() if n >= 5]
        # raw shared count favours long proteins; also rank by length-normalised share
        for key in (lambda t: t[1], lambda t: t[1] / math.sqrt(xk[t[0]])):
            for x, _ in sorted(hits, key=key, reverse=True)[:8]:
                xy_cand[y].add(x)
    sym2id = {g["symbol"]: gid for gid, g in mouse.items()}
    for y, x in XY_PAIRS:  # make sure the textbook pairs are tested
        if y in sym2id and x in sym2id:
            xy_cand[sym2id[y]].add(sym2id[x])
    pairs = [(y, x) for y, xs in xy_cand.items() for x in xs]
    with mp.Pool() as pool:
        xres = pool.map(compare, [(cds[y], cds[x], True) for y, x in pairs], chunksize=10)
    gam = {}
    for (y, x), r in zip(pairs, xres):
        if not r or r["cov"] < 0.5 or r["pid"] < 35:  # below this, end-gap-free hits are noise
            continue
        score = r["pid"] * min(r["cov"], 1)
        if y not in gam or score > gam[y]["score"]:
            gam[y] = {"x": x, "score": score, **r}

    # ---- 4. per-gene table
    genes = []
    for gid, g in mouse.items():
        loc = g["loc"] or {}
        chrom = loc.get("chr", "?")
        placed = chrom in ("X", "Y")
        orow = {}
        for sp in sp_name.values():
            o = best.get((gid, sp))
            if o:
                orow[sp] = o
        stratum = "mouse-specific"
        for name, sps in STRATA:
            if any(s in orow for s in sps):
                stratum = name
        go = g["go"]
        exp = sum(1 for t in go if t[3] in EXPERIMENTAL)
        fam = family.get(gid)
        rec = {
            "id": gid, "sym": g["symbol"], "name": g["name"],
            "chr": chrom if placed else ("Y" if "Y" in chrom else chrom),
            "placed": placed, "start": loc.get("start"), "end": loc.get("end"),
            "strand": loc.get("strand", ""),
            "aa": len(prot[gid]) if gid in prot else None,
            "stratum": stratum, "family": fam[0] if fam else None,
            "fam_size": fam[1] if fam else 1,
            "go_n": len(go), "go_exp": exp, "func": func_classes(go),
            "phast": None, "orth": orow,
        }
        if placed and loc.get("start"):
            m = phast_mean(pidx, chrom, loc["start"], loc["end"])
            rec["phast"] = None if m is None else round(m, 4)
        if gid in gam:
            x = gam[gid]
            rec["gametolog"] = {"x": mouse[x["x"]]["symbol"], "x_id": x["x"],
                                "x_start": (mouse[x["x"]]["loc"] or {}).get("start"),
                                "pid": x["pid"], "cov": x["cov"],
                                "dn": x.get("dn"), "ds": x.get("ds")}
        genes.append(rec)
    genes.sort(key=lambda r: (r["chr"], r["start"] or 0))

    # ---- 5. windows
    windows = []
    for chrom, L in CHROM_LEN.items():
        for w0 in range(0, L, WINDOW):
            w1 = min(w0 + WINDOW, L)
            gs = [r for r in genes if r["chr"] == chrom and r["placed"] and r["start"]
                  and w0 <= (r["start"] + r["end"]) // 2 < w1]
            row = {"chr": chrom, "start": w0, "end": w1, "n": len(gs),
                   "n_exp": sum(1 for r in gs if r["go_exp"] > 0),
                   "n_multi": sum(1 for r in gs if r["fam_size"] > 1),
                   "n_old": sum(1 for r in gs if r["stratum"] in ("Theria", "Mammalia", "Amniota",
                                                                   "Tetrapoda", "Vertebrata")),
                   "n_specific": sum(1 for r in gs if r["stratum"] in ("mouse-specific", "Murinae",
                                                                        "Rodentia")),
                   "phast": phast_mean(pidx, chrom, w0, w1)}
            for cls in FUNC_CLASSES:
                row["n_" + cls] = sum(1 for r in gs if cls in r["func"])
            for sp in ("rat", "human", "dog", "opossum", "platypus", "chicken"):
                v = [r["orth"][sp]["pid"] for r in gs if sp in r["orth"] and "pid" in r["orth"][sp]]
                row["pid_" + sp] = round(sum(v) / len(v), 1) if v else None
            v = [r["orth"]["rat"]["dn"] / r["orth"]["rat"]["ds"] for r in gs
                 if "rat" in r["orth"] and r["orth"]["rat"].get("ds") and r["orth"]["rat"].get("dn") is not None
                 and r["orth"]["rat"]["ds"] > 0.01]
            row["dnds_rat"] = round(sorted(v)[len(v) // 2], 3) if v else None
            if row["phast"] is not None:
                row["phast"] = round(row["phast"], 4)
            windows.append(row)

    # ---- 6. write tables
    species = [v[0] for v in PANEL.values()]
    with open(os.path.join(DATA, "genes.tsv"), "w") as f:
        cols = ["gene_id", "symbol", "chr", "start", "end", "aa", "phylostratum", "family",
                "family_size", "go_terms", "go_experimental", "functional_classes",
                "phastcons_gene_body", "x_gametolog", "xy_pid", "xy_dS",
                *[f"pid_{s}" for s in species], "dN_rat", "dS_rat", "dN_human", "dS_human"]
        f.write("\t".join(cols) + "\n")
        for r in genes:
            gm = r.get("gametolog") or {}
            vals = [r["id"], r["sym"], r["chr"], r["start"], r["end"], r["aa"], r["stratum"],
                    r["family"], r["fam_size"], r["go_n"], r["go_exp"], ",".join(r["func"]),
                    r["phast"], gm.get("x"), gm.get("pid"), gm.get("ds"),
                    *[r["orth"].get(s, {}).get("pid") for s in species],
                    r["orth"].get("rat", {}).get("dn"), r["orth"].get("rat", {}).get("ds"),
                    r["orth"].get("human", {}).get("dn"), r["orth"].get("human", {}).get("ds")]
            f.write("\t".join("" if v is None else str(v) for v in vals) + "\n")
    with open(os.path.join(DATA, "orthologs.tsv"), "w") as f:
        f.write("gene_id\tsymbol\tspecies\tclade\tortholog_gene_id\tortholog_symbol\tortholog_chr\t"
                "ortholog_start\tn_orthologs\tprotein_identity\tcoverage\tdN\tdS\tsource\n")
        for r in genes:
            for sp, o in r["orth"].items():
                vals = [r["id"], r["sym"], sp, sp_clade[sp], o["gid"], o["sym"], o["chr"],
                        o.get("start"), o["n"], o.get("pid"), o.get("cov"), o.get("dn"), o.get("ds"),
                        o.get("src", "NCBI")]
                f.write("\t".join("" if v is None else str(v) for v in vals) + "\n")
    with open(os.path.join(DATA, "gametologs.tsv"), "w") as f:
        f.write("y_gene_id\ty_symbol\ty_start\tfamily\tx_symbol\tx_start\tprotein_identity\t"
                "coverage\tdN\tdS\n")
        for r in genes:
            gm = r.get("gametolog")
            if gm:
                vals = [r["id"], r["sym"], r["start"], r["family"], gm["x"], gm["x_start"],
                        gm["pid"], gm["cov"], gm["dn"], gm["ds"]]
                f.write("\t".join("" if v is None else str(v) for v in vals) + "\n")
    with open(os.path.join(DATA, "windows.tsv"), "w") as f:
        cols = list(windows[0])
        f.write("\t".join(cols) + "\n")
        for w in windows:
            f.write("\t".join("" if w[c] is None else str(w[c]) for c in cols) + "\n")

    # ---- 7. compact JSON for the viewer
    slim_genes = []
    for r in genes:
        o = {sp: [v["chr"], v.get("start"), v.get("pid"), v["sym"], v.get("dn"), v.get("ds")]
             for sp, v in r["orth"].items()}
        g = {k: r[k] for k in ("id", "sym", "name", "chr", "placed", "start", "end", "aa", "stratum",
                               "family", "fam_size", "go_exp", "func", "phast")}
        g["o"] = o
        if "gametolog" in r:
            g["gam"] = r["gametolog"]
        slim_genes.append(g)
    with open(os.path.join(DATA, "xy_data.json"), "w") as f:
        json.dump({"species": [{"name": v[0], "clade": v[1], "taxid": t} for t, v in PANEL.items()],
                   "strata": [s for s, _ in STRATA], "chrom_len": CHROM_LEN, "window": WINDOW,
                   "genes": slim_genes, "windows": windows,
                   "phast_reduction": zoom["reduction"]},
                  f, separators=(",", ":"), ensure_ascii=False)
    print(f"{len(genes)} genes, {len(windows)} windows written", file=sys.stderr)


if __name__ == "__main__":
    main()
