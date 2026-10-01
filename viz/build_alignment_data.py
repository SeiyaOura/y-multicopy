#!/usr/bin/env python3
"""Build reference-anchored alignment data for the sequence viewer.

Each distinct CDS / protein sequence of a family is globally aligned to the
family reference (official symbol) and projected onto reference coordinates:
'.' = same as reference, letter = substitution, '-' = deletion,
insertions are listed separately.  Also computes pairwise identity between
all distinct sequences.
Output: viz/alignment_data.json and viz/index.html (viewer_template.html + data)
"""
import collections
import json
import os

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__))
SEQ = os.path.join(HERE, "..", "sequences")
FAMILIES = {
    "Ssty1": "spermiogenesis specific transcript on the Y 1",
    "Ssty2": "spermiogenesis specific transcript on the Y 2",
    "Sly": "Sycp3 like Y-linked",
    "Srsy": "serine-rich, secreted, Y-linked",
}


def aligner(kind):
    if kind == "aa":
        a = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                            open_gap_score=-11, extend_gap_score=-1)
    else:
        a = PairwiseAligner(mode="global", match_score=2, mismatch_score=-3,
                            open_gap_score=-6, extend_gap_score=-2)
    a.end_gap_score = 0  # truncated copies should not be penalised for missing ends
    return a


def project(ref, seq, al):
    """Return (line over ref coords, insertions, identity over aligned columns)."""
    a = al.align(ref, seq)[0]
    line = ["-"] * len(ref)
    ins = []
    match = cols = 0
    (rb, qb) = a.aligned
    if len(qb) and qb[0][0] > 0:  # N-terminal / 5' extension beyond the reference
        ins.append([0, seq[:qb[0][0]]])
    for i, ((r0, r1), (q0, q1)) in enumerate(zip(rb, qb)):
        for k in range(r1 - r0):
            c = seq[q0 + k]
            line[r0 + k] = "." if c == ref[r0 + k] else c
            match += c == ref[r0 + k]
            cols += 1
        if i + 1 < len(rb):
            nr0, nq0 = rb[i + 1][0], qb[i + 1][0]
            cols += (nr0 - r1) + (nq0 - q1)  # internal gaps count against identity
            if nq0 > q1:
                ins.append([int(r1), seq[q1:nq0]])
    if len(qb) and qb[-1][1] < len(seq):  # C-terminal / 3' extension
        ins.append([len(ref), seq[qb[-1][1]:]])
    return "".join(line), ins, float(match / cols) if cols else 0.0


def identity(a_seq, b_seq, al):
    _, _, pid = project(a_seq, b_seq, al)
    return round(pid, 4)


def family_data(fam):
    meta = {}
    for row in __import__("csv").DictReader(open(f"{SEQ}/{fam}.genes.tsv"), delimiter="\t"):
        meta[row["protein"]] = row
    out = {}
    for kind, path in (("aa", f"{SEQ}/{fam}.protein.faa"), ("nt", f"{SEQ}/{fam}.cds.fna")):
        groups = collections.OrderedDict()
        for r in SeqIO.parse(path, "fasta"):
            prot = r.description.split()[1]
            groups.setdefault(str(r.seq), []).append(prot)
        seqs = list(groups)
        ref = seqs[0]  # files list the official-symbol reference first
        order = [0] + sorted(range(1, len(seqs)), key=lambda i: -len(groups[seqs[i]]))
        al = aligner(kind)
        rows = []
        for n, i in enumerate(order):
            s, prots = seqs[i], groups[seqs[i]]
            line, ins, pid = project(ref, s, al)
            genes = sorted({meta[p]["symbol"] for p in prots})
            rows.append(dict(name="ref" if n == 0 else f"v{n}", copies=len(prots), length=len(s),
                             identity=round(pid, 4), genes=genes, proteins=prots,
                             line=line, ins=ins, seq=s))
        ordered = [seqs[i] for i in order]
        matrix = [[1.0 if a == b else identity(ordered[a], ordered[b], al) if a < b else None
                   for b in range(len(ordered))] for a in range(len(ordered))]
        for a in range(len(ordered)):
            for b in range(a):
                matrix[a][b] = matrix[b][a]
        out[kind] = dict(ref=ref, rows=rows, matrix=matrix)
    return out


def main():
    data = {"families": {}}
    for fam, desc in FAMILIES.items():
        d = family_data(fam)
        d["description"] = desc
        d["transcripts"] = sum(r["copies"] for r in d["aa"]["rows"])
        d["genes"] = len({g for r in d["aa"]["rows"] for g in r["genes"]})
        data["families"][fam] = d
        print(fam, len(d["aa"]["rows"]), "protein variants,", len(d["nt"]["rows"]), "CDS variants")
    # Ssty1 vs Ssty2 paralog comparison, reference proteins only
    s1 = data["families"]["Ssty1"]["aa"]["ref"]
    s2 = data["families"]["Ssty2"]["aa"]["ref"]
    a = aligner("aa").align(s1, s2)[0]
    data["paralog"] = dict(a=str(a[0]), b=str(a[1]), a_name="Ssty1", b_name="Ssty2",
                           identity=identity(s1, s2, aligner("aa")))
    blob = json.dumps(data, separators=(",", ":"))
    with open(os.path.join(HERE, "alignment_data.json"), "w") as f:
        f.write(blob)
    # inline the data into the standalone viewer page
    with open(os.path.join(HERE, "viewer_template.html")) as f:
        page = f.read().replace("/*DATA*/null", blob)
    with open(os.path.join(HERE, "index.html"), "w") as f:
        f.write(page)


if __name__ == "__main__":
    main()
