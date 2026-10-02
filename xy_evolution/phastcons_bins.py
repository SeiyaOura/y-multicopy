#!/usr/bin/env python3
"""Summarise UCSC mm39 phastCons35way (35 vertebrates) over chrX / chrY
without downloading the 7 GB bigWig: read one zoom level through HTTP range
requests (pure python, no libcurl needed).

Output: data/phastcons_zoom.json = {chrom: [[start, end, validCount, sum], ...]}
These records are aggregated per window / per gene in build_xy_data.py.
"""
import json
import os
import struct
import sys
import urllib.request
import zlib

URL = "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/phastCons35way/mm39.phastCons35way.bw"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data", "phastcons_zoom.json")
CHROMS = ("chrX", "chrY")
MAX_REDUCTION = 5000  # finest zoom level with reduction <= this (bp)


def fetch(off, size):
    req = urllib.request.Request(URL, headers={"Range": f"bytes={off}-{off + size - 1}"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read()
        except OSError:
            if attempt == 3:
                raise


def main():
    hdr = fetch(0, 64 + 24 * 16)
    magic, ver, nzoom, ctree, data_off, idx_off = struct.unpack_from("<IHHQQQ", hdr, 0)
    assert magic == 0x888FFC26, "not a bigWig"
    uncompress = struct.unpack_from("<I", hdr, 52)[0]
    zooms = [struct.unpack_from("<IIQQ", hdr, 64 + 24 * i) for i in range(nzoom)]
    level = max((z for z in zooms if z[0] <= MAX_REDUCTION), key=lambda z: z[0])
    print(f"zoom levels {[z[0] for z in zooms]}; using {level[0]} bp", file=sys.stderr)

    # chromosome B+ tree
    buf = fetch(ctree, 32)
    _, bsize, ksize, vsize, nitems = struct.unpack_from("<IIIIQ", buf, 0)
    chrom_ids = {}

    def walk_chrom(off):
        head = fetch(off, 4)
        leaf, _, cnt = struct.unpack("<BBH", head)
        item = ksize + (8 if not leaf else vsize)
        body = fetch(off + 4, cnt * item)
        for i in range(cnt):
            key = body[i * item:i * item + ksize].rstrip(b"\0").decode()
            if leaf:
                cid, csize = struct.unpack_from("<II", body, i * item + ksize)
                chrom_ids[key] = (cid, csize)
            else:
                walk_chrom(struct.unpack_from("<Q", body, i * item + ksize)[0])

    walk_chrom(ctree + 32)
    want = {chrom_ids[c][0]: c for c in CHROMS}

    # R-tree of the chosen zoom level -> leaf blocks overlapping chrX/chrY
    blocks = []

    def walk_r(off):
        head = fetch(off, 4)
        leaf, _, cnt = struct.unpack("<BBH", head)
        item = 32 if leaf else 24
        body = fetch(off + 4, cnt * item)
        for i in range(cnt):
            sc, sb, ec, eb, child = struct.unpack_from("<IIIIQ", body, i * item)
            if not any(sc <= c <= ec for c in want):
                continue
            if leaf:
                blocks.append((child, struct.unpack_from("<Q", body, i * item + 24)[0]))
            else:
                walk_r(child)

    walk_r(level[3] + 48)
    print(f"{len(blocks)} zoom blocks", file=sys.stderr)

    out = {c: [] for c in CHROMS}
    blocks.sort()
    # merge adjacent blocks into larger range requests
    i = 0
    while i < len(blocks):
        start, end = blocks[i][0], blocks[i][0] + blocks[i][1]
        j = i + 1
        while j < len(blocks) and blocks[j][0] == end and end - start < 8_000_000:
            end += blocks[j][1]
            j += 1
        raw = fetch(start, end - start)
        for off, size in blocks[i:j]:
            data = raw[off - start:off - start + size]
            if uncompress:
                data = zlib.decompress(data)
            for k in range(0, len(data), 32):
                cid, s, e, n, mn, mx, sm, sq = struct.unpack_from("<IIIIffff", data, k)
                if cid in want:
                    out[want[cid]].append([s, e, n, round(sm, 3)])
        i = j
    for c in CHROMS:
        out[c].sort()
        print(f"{c}: {len(out[c])} records", file=sys.stderr)
    with open(OUT, "w") as f:
        json.dump({"reduction": level[0], "source": URL, **out}, f, separators=(",", ":"))


if __name__ == "__main__":
    main()
