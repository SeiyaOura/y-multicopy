#!/usr/bin/env python3
"""Fetch mouse chrX / chrY protein-coding genes and their orthologs across a
vertebrate panel from NCBI (Gene + Datasets v2 API).

Steps
  1. list every live protein-coding mouse gene on X and Y (Entrez esearch),
  2. gene reports (GRCm39 position, GO annotations with evidence codes),
  3. NCBI Orthologs for each mouse gene, restricted to the species panel
     (chromosome + position of each ortholog),
  4. CDS of every mouse gene and every ortholog (longest CDS per gene kept),
  5. CDS of the X-linked gametologs of mouse Y genes (for X-Y divergence).

Output: data/cache/*.json, data/cache/cds.fna.gz  (re-runs reuse the cache)
Usage:  python3 fetch_xy_data.py          (optional env NCBI_API_KEY)
"""
import gzip
import io
import json
import os
import sys
import time
import urllib.parse
import urllib.error
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "data", "cache")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/"
MOUSE = 10090

# taxid: (short name, clade) - ordered from closest to most distant
PANEL = {
    10116: ("rat", "Rodentia"),
    10036: ("hamster", "Rodentia"),
    10141: ("guinea_pig", "Rodentia"),
    9606: ("human", "Primates"),
    9598: ("chimp", "Primates"),
    9544: ("macaque", "Primates"),
    9483: ("marmoset", "Primates"),
    9615: ("dog", "Laurasiatheria"),
    9913: ("cow", "Laurasiatheria"),
    9796: ("horse", "Laurasiatheria"),
    9785: ("elephant", "Afrotheria"),
    13616: ("opossum", "Marsupialia"),
    9258: ("platypus", "Monotremata"),
    9031: ("chicken", "Aves"),
    28377: ("anole", "Squamata"),
    8364: ("xenopus", "Amphibia"),
    7955: ("zebrafish", "Teleostei"),
}

ENSEMBL = "https://rest.ensembl.org/"
ENS_SPECIES = {
    10116: "rattus_norvegicus", 10036: "mesocricetus_auratus", 10141: "cavia_porcellus",
    9606: "homo_sapiens", 9598: "pan_troglodytes", 9544: "macaca_mulatta", 9483: "callithrix_jacchus",
    9615: "canis_lupus_familiaris", 9913: "bos_taurus", 9796: "equus_caballus",
    9785: "loxodonta_africana", 13616: "monodelphis_domestica", 9258: "ornithorhynchus_anatinus",
    9031: "gallus_gallus", 28377: "anolis_carolinensis", 8364: "xenopus_tropicalis", 7955: "danio_rerio",
}

# X-linked gametologs (or X homologs) of mouse Y genes, by NCBI symbol
XY_PAIRS = [
    ("Sry", "Sox3"), ("Rbmy", "Rbmx"), ("Zfy1", "Zfx"), ("Zfy2", "Zfx"),
    ("Uba1y", "Uba1"), ("Kdm5d", "Kdm5c"), ("Eif2s3y", "Eif2s3x"),
    ("Uty", "Kdm6a"), ("Ddx3y", "Ddx3x"), ("Usp9y", "Usp9x"),
    ("H2al2y", "H2al2a"), ("Sly", "Slx"), ("Ssty1", "Sstx"), ("Ssty2", "Sstx"),
    ("Srsy", "Srsx"),
]


def _open(url, data=None, headers=None, timeout=600):
    req = urllib.request.Request(url, data=data, headers=headers or {})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                body = r.read()
            time.sleep(0.12 if os.environ.get("NCBI_API_KEY") else 0.25)
            return body
        except OSError as e:
            if attempt == 4:
                raise
            print(f"  retry {url[:80]} ({e})", file=sys.stderr)
            time.sleep(2 ** (attempt + 1))


def datasets(path, payload=None):
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    if os.environ.get("NCBI_API_KEY"):
        headers["api-key"] = os.environ["NCBI_API_KEY"]
    data = json.dumps(payload).encode() if payload is not None else None
    return _open(DATASETS + path, data, headers)


def cached(name, fn):
    path = os.path.join(CACHE, name)
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    obj = fn()
    with open(path, "w") as f:
        json.dump(obj, f)
    return obj


def esearch_ids(term):
    params = {"db": "gene", "term": term, "retmax": 100000, "retmode": "json"}
    if os.environ.get("NCBI_API_KEY"):
        params["api_key"] = os.environ["NCBI_API_KEY"]
    body = _open(EUTILS + "esearch.fcgi", urllib.parse.urlencode(params).encode())
    return json.loads(body)["esearchresult"]["idlist"]


def gene_reports(ids):
    out = {}
    for i in range(0, len(ids), 500):
        payload = {"gene_ids": [int(x) for x in ids[i:i + 500]], "page_size": 1000}
        while True:
            res = json.loads(datasets("gene", payload))
            for r in res.get("reports", []):
                out[r["gene"]["gene_id"]] = slim(r["gene"])
            if not res.get("next_page_token"):
                break
            payload["page_token"] = res["next_page_token"]
    return out


def slim(g):
    """Keep only the fields used downstream."""
    loc = None
    for ann in g.get("annotations", []):
        for gl in ann.get("genomic_locations", []):
            loc = {"chr": gl["sequence_name"], "start": int(gl["genomic_range"]["begin"]),
                   "end": int(gl["genomic_range"]["end"]),
                   "strand": gl["genomic_range"].get("orientation", "")}
            break
        if loc:
            break
    go = []
    for kind, terms in (g.get("gene_ontology") or {}).items():
        if isinstance(terms, list):
            for t in terms:
                go.append([kind[0].upper(), t["go_id"], t["name"], t.get("evidence_code", "")])
    return {"gene_id": g["gene_id"], "symbol": g["symbol"], "name": g.get("description", ""),
            "tax_id": int(g["tax_id"]), "chromosomes": g.get("chromosomes", []), "loc": loc,
            "ens": g.get("ensembl_gene_ids", []), "go": go}


def orthologs(gene_id):
    q = "&".join(f"taxon_filter={t}" for t in PANEL)
    out, token = [], ""
    while True:
        res = json.loads(datasets(f"gene/id/{gene_id}/orthologs?page_size=1000&{q}{token}"))
        out += [slim(r["gene"]) for r in res.get("reports", []) if int(r["gene"]["tax_id"]) in PANEL]
        if not res.get("next_page_token"):
            return out
        token = "&page_token=" + urllib.parse.quote(res["next_page_token"])


def ensembl(path, payload=None):
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    data = json.dumps(payload).encode() if payload is not None else None
    try:
        return json.loads(_open(ENSEMBL + path, data, headers))
    except urllib.error.HTTPError as e:
        if e.code in (400, 404):  # unknown stable id
            return None
        raise


def ensembl_orthologs(ens_id):
    """Ensembl Compara orthologs of one mouse gene in the panel species."""
    q = ";".join(f"target_species={s}" for s in ENS_SPECIES.values())
    res = ensembl(f"homology/id/mus_musculus/{ens_id}?type=orthologues;sequence=none;{q}")
    out = []
    for d in (res or {}).get("data", []):
        for h in d["homologies"]:
            out.append({"id": h["target"]["id"], "species": h["target"]["species"],
                        "type": h["type"], "pid_target": h["target"]["perc_id"],
                        "pid_source": h["source"]["perc_id"], "dn_ds": h.get("dn_ds")})
    return out


def ensembl_locations(ids):
    loc = {}
    ids = sorted(ids)
    for i in range(0, len(ids), 1000):
        res = ensembl("lookup/id", {"ids": ids[i:i + 1000]}) or {}
        for k, v in res.items():
            if v:
                loc[k] = {"chr": v.get("seq_region_name"), "start": v.get("start"),
                          "sym": v.get("display_name") or k}
    return loc


def longest_cds(ids):
    """Download CDS FASTA for gene ids; return {gene_id: longest CDS (str)}."""
    best = {}
    for i in range(0, len(ids), 400):
        chunk = [int(x) for x in ids[i:i + 400]]
        z = zipfile.ZipFile(io.BytesIO(datasets(
            "gene/download", {"gene_ids": chunk, "include_annotation_type": ["FASTA_CDS"]})))
        names = [n for n in z.namelist() if n.endswith("cds.fna")]
        if not names:
            continue
        gid, seq = None, []
        for line in io.TextIOWrapper(z.open(names[0])):
            if line.startswith(">"):
                _keep(best, gid, seq)
                gid = line.split("[GeneID=")[1].split("]")[0] if "[GeneID=" in line else None
                seq = []
            else:
                seq.append(line.strip())
        _keep(best, gid, seq)
        print(f"  CDS {min(i + 400, len(ids))}/{len(ids)}", file=sys.stderr)
    return best


def _keep(best, gid, seq):
    if gid is None:
        return
    s = "".join(seq).upper()
    if len(s) % 3 == 0 and len(s) > len(best.get(gid, "")):
        best[gid] = s


def main():
    os.makedirs(CACHE, exist_ok=True)
    ids = {}
    for c in ("X", "Y"):
        ids[c] = cached(f"mouse_{c}_ids.json", lambda c=c: esearch_ids(
            f"{MOUSE}[taxid] AND {c}[chr] AND alive[prop] AND genetype protein coding[prop]"))
        print(f"mouse chr{c}: {len(ids[c])} protein-coding genes", file=sys.stderr)
    all_ids = ids["X"] + ids["Y"]
    mouse = cached("mouse_genes.json", lambda: gene_reports(all_ids))

    # X gametologs of Y genes that may lie elsewhere (all expected on X, but be safe)
    sym2id = {g["symbol"]: gid for gid, g in mouse.items()}
    missing = [x for _, x in XY_PAIRS if x not in sym2id]
    if missing:
        print("gametolog symbols not found on X:", missing, file=sys.stderr)

    orth_path = os.path.join(CACHE, "orthologs.json")
    orth = json.load(open(orth_path)) if os.path.exists(orth_path) else {}
    todo = [g for g in all_ids if g not in orth]
    with ThreadPoolExecutor(4) as pool:  # ~4 requests in flight stays under 5 req/s
        for n, (gid, res) in enumerate(zip(todo, pool.map(orthologs, todo)), 1):
            orth[gid] = res
            if n % 50 == 0 or n == len(todo):
                print(f"  orthologs {n}/{len(todo)}", file=sys.stderr)
                with open(orth_path, "w") as f:
                    json.dump(orth, f)

    ens_path = os.path.join(CACHE, "ensembl_orthologs.json")
    ens = json.load(open(ens_path)) if os.path.exists(ens_path) else {}
    todo = [(gid, e) for gid in all_ids if gid in mouse for e in mouse[gid]["ens"] if e not in ens]
    with ThreadPoolExecutor(4) as pool:
        for n, ((gid, e), res) in enumerate(zip(todo, pool.map(lambda t: ensembl_orthologs(t[1]), todo)), 1):
            ens[e] = res
            if n % 100 == 0 or n == len(todo):
                print(f"  Ensembl orthologs {n}/{len(todo)}", file=sys.stderr)
                with open(ens_path, "w") as f:
                    json.dump(ens, f)
    cached("ensembl_locations.json", lambda: ensembl_locations(
        {h["id"] for v in ens.values() for h in v}))

    cds_path = os.path.join(CACHE, "cds.fna.gz")
    want = set(all_ids) | {o["gene_id"] for v in orth.values() for o in v}
    have = {}
    if os.path.exists(cds_path):
        with gzip.open(cds_path, "rt") as f:
            for line in f:
                if line.startswith(">"):
                    gid = line[1:].strip()
                else:
                    have[gid] = line.strip()
    new = longest_cds(sorted(want - set(have)))
    have.update(new)
    with gzip.open(cds_path, "wt") as f:
        for gid, s in have.items():
            f.write(f">{gid}\n{s}\n")
    print(f"CDS for {len(have)} / {len(want)} genes", file=sys.stderr)


if __name__ == "__main__":
    main()
