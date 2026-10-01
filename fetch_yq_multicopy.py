#!/usr/bin/env python3
"""Fetch every RefSeq CDS (DNA) and protein sequence for the mouse Yq
multicopy gene families Ssty1, Ssty2, Sly and Srsy from NCBI.

Many copies on Yq are annotated under Gm/LOC symbols, so genes are matched
both by official symbol and by full gene name / description.

Usage: python3 fetch_yq_multicopy.py [outdir]
Output: <outdir>/<family>.cds.fna, <family>.protein.faa, <family>.genes.tsv
Requires network access to eutils.ncbi.nlm.nih.gov. Standard library only.
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
ORG = "10090[taxid] AND Y[chr] AND alive[prop]"

FAMILIES = {
    "Ssty1": '(Ssty1[sym] OR "spermiogenesis specific transcript on the Y 1"[All Fields])',
    "Ssty2": '(Ssty2[sym] OR "spermiogenesis specific transcript on the Y 2"[All Fields])',
    "Sly": '(Sly[sym] OR "Sycp3 like Y-linked"[All Fields])',
    "Srsy": '(Srsy[sym] OR "serine-rich, secreted, Y-linked"[All Fields])',
}


def call(tool, **params):
    params.setdefault("tool", "yq_multicopy")
    if os.environ.get("NCBI_API_KEY"):
        params["api_key"] = os.environ["NCBI_API_KEY"]
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(4):
        try:
            with urllib.request.urlopen(EUTILS + tool, data=data, timeout=120) as r:
                body = r.read().decode()
            time.sleep(0.4)  # stay under 3 requests/s without an API key
            return body
        except OSError:
            if attempt == 3:
                raise
            time.sleep(2 ** (attempt + 1))


def gene_ids(term):
    res = json.loads(call("esearch.fcgi", db="gene", term=f"{term} AND {ORG}",
                          retmax=1000, retmode="json"))
    return res["esearchresult"]["idlist"]


def gene_table(ids):
    res = json.loads(call("esummary.fcgi", db="gene", id=",".join(ids), retmode="json"))
    rows = []
    for gid in res["result"]["uids"]:
        g = res["result"][gid]
        loc = (g.get("genomicinfo") or [{}])[0]
        rows.append([gid, g["name"], g["description"], g.get("chromosome", ""),
                     str(loc.get("chrstart", "")), str(loc.get("chrstop", ""))])
    return rows


def refseq_mrnas(ids):
    res = json.loads(call("elink.fcgi", dbfrom="gene", db="nuccore", id=",".join(ids),
                          linkname="gene_nuccore_refseqrna", retmode="json"))
    out = []
    for ls in res["linksets"]:
        for db in ls.get("linksetdbs", []):
            out += db["links"]
    return sorted(set(out))


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "sequences"
    os.makedirs(outdir, exist_ok=True)
    for fam, term in FAMILIES.items():
        ids = gene_ids(term)
        print(f"{fam}: {len(ids)} genes", file=sys.stderr)
        if not ids:
            continue
        with open(f"{outdir}/{fam}.genes.tsv", "w") as f:
            f.write("gene_id\tsymbol\tdescription\tchr\tstart\tstop\n")
            for row in gene_table(ids):
                f.write("\t".join(row) + "\n")
        mrnas = refseq_mrnas(ids)
        print(f"{fam}: {len(mrnas)} RefSeq transcripts", file=sys.stderr)
        if not mrnas:
            continue
        for rettype, ext in (("fasta_cds_na", "cds.fna"), ("fasta_cds_aa", "protein.faa")):
            fasta = call("efetch.fcgi", db="nuccore", id=",".join(mrnas),
                         rettype=rettype, retmode="text")
            with open(f"{outdir}/{fam}.{ext}", "w") as f:
                f.write(fasta)


if __name__ == "__main__":
    main()
