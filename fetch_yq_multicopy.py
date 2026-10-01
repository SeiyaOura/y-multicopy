#!/usr/bin/env python3
"""Collect every RefSeq CDS (DNA) and protein sequence of the mouse Yq
multicopy families Ssty1, Ssty2, Sly and Srsy.

Most Yq copies are annotated under Gm/LOC symbols, so symbol search is not
enough.  Instead:
  1. fetch the reference protein of each family (official symbol),
  2. fetch the CDS/protein of every RefSeq transcript of every live mouse
     chrY gene,
  3. assign each protein to the family whose reference gives the best local
     alignment score (BLOSUM62), keeping hits with score >= MIN_SCORE.

Usage: python3 fetch_yq_multicopy.py [outdir]     (needs: pip install biopython)
Output per family: <family>.cds.fna, .protein.faa, .protein.unique.faa, .genes.tsv
"""
import collections
import io
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
FAMILIES = ["Ssty1", "Ssty2", "Sly", "Srsy"]
MIN_SCORE = 300


def call(tool, **params):
    params.setdefault("tool", "yq_multicopy")
    if os.environ.get("NCBI_API_KEY"):
        params["api_key"] = os.environ["NCBI_API_KEY"]
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(4):
        try:
            with urllib.request.urlopen(EUTILS + tool, data=data, timeout=300) as r:
                body = r.read().decode()
            time.sleep(0.4)  # stay under 3 requests/s without an API key
            return body
        except OSError:
            if attempt == 3:
                raise
            time.sleep(2 ** (attempt + 1))


def esearch(term):
    res = json.loads(call("esearch.fcgi", db="gene", term=term, retmax=10000, retmode="json"))
    return res["esearchresult"]["idlist"]


def gene_info(ids):
    info = {}
    for i in range(0, len(ids), 400):
        res = json.loads(call("esummary.fcgi", db="gene", id=",".join(ids[i:i + 400]),
                              retmode="json"))["result"]
        for gid in res["uids"]:
            g = res[gid]
            loc = (g.get("genomicinfo") or [{}])[0]
            info[g["name"]] = (gid, str(loc.get("chrstart", "")))
    return info


def refseq_rnas(ids):
    acc = set()
    for i in range(0, len(ids), 200):
        res = json.loads(call("elink.fcgi", dbfrom="gene", db="nuccore", id=",".join(ids[i:i + 200]),
                              linkname="gene_nuccore_refseqrna", retmode="json"))
        for ls in res["linksets"]:
            for db in ls.get("linksetdbs", []):
                acc.update(db["links"])
    return sorted(acc)


def efetch(acc, rettype):
    text = "".join(call("efetch.fcgi", db="nuccore", id=",".join(acc[i:i + 200]),
                        rettype=rettype, retmode="text") for i in range(0, len(acc), 200))
    return list(SeqIO.parse(io.StringIO(text), "fasta"))


def tag(rec, key):
    m = re.search(rf"\[{key}=([^\]]+)\]", rec.description)
    return m.group(1) if m else ""


def wrap(s):
    return "\n".join(s[i:i + 70] for i in range(0, len(s), 70))


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "sequences"
    os.makedirs(outdir, exist_ok=True)

    refs = {}
    for fam in FAMILIES:
        ids = esearch(f"{fam}[sym] AND 10090[taxid]")
        aa = efetch(refseq_rnas(ids), "fasta_cds_aa")
        refs[fam] = str(next(r for r in aa if tag(r, "gene") == fam).seq)

    ids = esearch("10090[taxid] AND Y[chr] AND alive[prop]")
    info = gene_info(ids)
    acc = refseq_rnas(ids)
    print(f"chrY: {len(ids)} genes, {len(acc)} RefSeq transcripts", file=sys.stderr)
    aa = efetch(acc, "fasta_cds_aa")
    na = {tag(r, "protein_id"): r for r in efetch(acc, "fasta_cds_na")}

    aligner = PairwiseAligner(mode="local", open_gap_score=-11, extend_gap_score=-1,
                              substitution_matrix=substitution_matrices.load("BLOSUM62"))
    hits = collections.defaultdict(list)
    for r in aa:
        seq = str(r.seq).rstrip("*")
        best = None
        for fam, ref in refs.items():
            a = aligner.align(seq, ref)[0]
            alen = sum(e - s for s, e in a.aligned[0])
            ident = sum(seq[s1 + k] == ref[s2 + k]
                        for (s1, e1), (s2, _) in zip(*a.aligned) for k in range(e1 - s1))
            cand = (a.score, fam, ident / alen if alen else 0, alen / len(ref))
            best = max(best, cand) if best else cand
        if best[0] < MIN_SCORE:
            continue
        gene, prot = tag(r, "gene"), tag(r, "protein_id")
        gid, start = info.get(gene, ("", ""))
        tx = r.id.split("|")[1].split("_prot_")[0]
        hits[best[1]].append(dict(gene=gene, gid=gid, start=start, tx=tx, prot=prot, aa=seq,
                                  na=str(na[prot].seq), pid=best[2], cov=best[3]))

    for fam in FAMILIES:
        rows = sorted(hits[fam], key=lambda x: (x["gene"] != fam, int(x["start"] or 0), x["tx"]))
        uniq = collections.OrderedDict()
        for x in rows:
            uniq.setdefault(x["aa"], []).append(x["prot"])
        with open(f"{outdir}/{fam}.cds.fna", "w") as fna, \
                open(f"{outdir}/{fam}.protein.faa", "w") as faa, \
                open(f"{outdir}/{fam}.genes.tsv", "w") as tsv:
            tsv.write("symbol\tgene_id\tchrY_start(GRCm39)\ttranscript\tprotein\tcds_nt"
                      "\tprotein_aa\tidentity_to_ref\tcoverage_of_ref\n")
            for x in rows:
                h = f"{x['tx']} {x['prot']} gene={x['gene']} GeneID={x['gid']} chrY:{x['start']}"
                fna.write(f">{h}\n{wrap(x['na'])}\n")
                faa.write(f">{h}\n{wrap(x['aa'])}\n")
                tsv.write(f"{x['gene']}\t{x['gid']}\t{x['start']}\t{x['tx']}\t{x['prot']}"
                          f"\t{len(x['na'])}\t{len(x['aa'])}\t{x['pid']:.3f}\t{x['cov']:.3f}\n")
        with open(f"{outdir}/{fam}.protein.unique.faa", "w") as f:
            for i, (s, prots) in enumerate(uniq.items(), 1):
                f.write(f">{fam}_variant{i} len={len(s)} n={len(prots)} {','.join(prots)}\n{wrap(s)}\n")
        print(f"{fam}: {len(rows)} transcripts, {len({x['gene'] for x in rows})} genes, "
              f"{len(uniq)} distinct proteins", file=sys.stderr)


if __name__ == "__main__":
    main()
