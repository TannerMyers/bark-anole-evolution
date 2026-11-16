#!/usr/bin/env python3
"""
================================================================================
DISTICHUS CHROMOSOMAL DISCORDANCE ANALYSIS PIPELINE
================================================================================

Comprehensive analysis script to compute phylogenetic discordance metrics,
support statistics, clade monophyly, stray destinations, and backbone RF
distances across four genomic partitions:
  1. Autosomes (90% completeness ddRAD)
  2. Y chromosome (50% completeness ddRAD)
  3. Ancestral X chromosome (X1, 50% completeness ddRAD)
  4. Neo-X chromosome (X2, 50% completeness ddRAD)

Also triggers reproduction of publication-grade Figure 2 style tanglegrams:
  - discordance_tanglegram_Y_auto.pdf  (Autosomes vs. Y)
  - discordance_tanglegram_x1_x2.pdf   (X1 vs. X2)

Usage:
  python3 calculate_discordance_stats_and_figures.py
  python3 calculate_discordance_stats_and_figures.py --skip-figures
  python3 calculate_discordance_stats_and_figures.py --outdir ./results
================================================================================
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

# ------------------------------------------------------------------------------
# 1. FILE RESOLUTION HELPERS
# ------------------------------------------------------------------------------
def find_file(candidates):
    for path in candidates:
        if os.path.exists(path):
            return os.path.abspath(path)
    return None

# ------------------------------------------------------------------------------
# 2. TREE PARSING & ANNOTATION (STANDALONE / NO EXTERNAL LIBRARIES)
# ------------------------------------------------------------------------------
class Node:
    __slots__ = ("label", "children", "brlen", "parent")

    def __init__(self, label=None, brlen=None):
        self.label = label
        self.children = []
        self.brlen = brlen
        self.parent = None


def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = [0]

    def node():
        n = Node()
        if s[pos[0]] == "(":
            pos[0] += 1
            while True:
                c = node()
                n.children.append(c)
                c.parent = n
                if s[pos[0]] == ",":
                    pos[0] += 1
                elif s[pos[0]] == ")":
                    pos[0] += 1
                    break
        m = re.match(r"[^():,;]*", s[pos[0]:])
        if m and m.group(0):
            n.label = m.group(0)
            pos[0] += len(m.group(0))
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            m2 = re.match(r"[0-9.eE+-]+", s[pos[0]:])
            if m2:
                n.brlen = float(m2.group(0))
                pos[0] += len(m2.group(0))
        return n

    t = node()
    assert pos[0] == len(s), f"Leftover characters at position {pos[0]}: {s[pos[0]:]}"
    return t


def norm_leaf(lbl):
    return lbl.split(".")[0]


def get_leaves(t):
    if not t.children:
        return [t.label]
    out = []
    for c in t.children:
        out += get_leaves(c)
    return out


def parse_support(label):
    if not label:
        return None
    m = re.match(r"^([0-9.]+)/([0-9.]+)$", label.strip())
    if m:
        return (float(m.group(1)), float(m.group(2)))
    m2 = re.match(r"^[0-9.]+$", label.strip())
    if m2:
        return (float(label), None)
    return None


def annotate(t, tip_to_tax):
    info = {}
    ctr = [0]

    def dfs(n, parent):
        nid = ctr[0]
        ctr[0] += 1
        if not n.children:
            leaf = norm_leaf(n.label)
            s = frozenset([leaf])
            ts = frozenset([tip_to_tax.get(leaf, "?")])
        else:
            res = [dfs(c, nid) for c in n.children]
            s = frozenset().union(*(r[0] for r in res))
            ts = frozenset().union(*(r[1] for r in res))
        info[nid] = {
            "label": n.label,
            "leafset": s,
            "taxa": ts,
            "parent": parent,
            "brlen": n.brlen,
        }
        return (s, ts)

    dfs(t, None)
    return info


def node_by_id(t):
    m = {}
    ctr = [0]

    def dfs(n):
        nid = ctr[0]
        ctr[0] += 1
        m[nid] = n
        for c in n.children:
            dfs(c)

    dfs(t)
    return m


def reroot_edge(t, x):
    p = x.parent
    if p is None:
        return t
    path = []
    cur = p
    while cur is not None:
        path.append(cur)
        cur = cur.parent
    r = Node()
    p.children.remove(x)
    x.parent = r
    p.parent = r
    r.children = [x, p]
    prev = p
    for i in range(1, len(path)):
        cur = path[i]
        cur.children.remove(prev)
        prev.children.append(cur)
        cur.parent = prev
        prev = cur
    if path[-1].brlen is None:
        path[-1].brlen = 0.0
    return r


def smallest_clade_node(t, seed_norms, tip_to_tax):
    info = annotate(t, tip_to_tax)
    nodes = node_by_id(t)
    idx = {}
    for nid, rec in info.items():
        if len(rec["leafset"]) == 1 and rec["label"]:
            l = next(iter(rec["leafset"]))
            if l not in idx:
                idx[l] = nid
    target = set(seed_norms)
    missing = [s for s in target if s not in idx]
    if missing:
        raise ValueError(f"Seed tips missing from tree: {missing}")
    cur = idx[seed_norms[0]]
    while not (target <= info[cur]["leafset"]):
        p = info[cur]["parent"]
        if p is None:
            raise ValueError("No clade contains all seeds")
        cur = p
    return cur, info, nodes

# ------------------------------------------------------------------------------
# 3. SPLITS, SUPPORT, & ROBINSON-FOULDS
# ------------------------------------------------------------------------------
def splits_with_support(info, all_tips, tip_to_tax):
    out = {}
    for nid, d in info.items():
        s = d["leafset"]
        if len(s) in (0, 1, len(all_tips)):
            continue
        a, b = s, frozenset(all_tips - s)
        small = a if len(a) <= len(b) else b
        if len(small) == 1:
            continue
        sup = parse_support(d["label"])
        taxa = frozenset(tip_to_tax.get(x, "?") for x in small)
        out[small] = {"support": sup, "taxa": taxa, "n": len(small)}
    return out


def support_stats(splits):
    vals_a = [v["support"][0] for v in splits.values() if v["support"]]
    vals_b = [v["support"][1] for v in splits.values() if v["support"] and v["support"][1] is not None]
    n = len(splits)
    return {
        "n_splits": n,
        "n_labelled": len(vals_a),
        "shalrt_mean": round(sum(vals_a) / len(vals_a), 1) if vals_a else None,
        "shalrt_median": sorted(vals_a)[len(vals_a) // 2] if vals_a else None,
        "ufboot_mean": round(sum(vals_b) / len(vals_b), 1) if vals_b else None,
        "ufboot_median": sorted(vals_b)[len(vals_b) // 2] if vals_b else None,
        "pct_shalrt_ge95": round(100 * sum(1 for v in vals_a if v >= 95) / len(vals_a), 1) if vals_a else None,
        "pct_ufboot_ge95": round(100 * sum(1 for v in vals_b if v >= 95) / len(vals_b), 1) if vals_b else None,
        "pct_both100": round(100 * sum(1 for v in vals_a if v == 100) / len(vals_a), 1) if vals_a else None,
    }


def compute_rf(splits_a, splits_b, n_tips=182):
    set_a = set(splits_a.keys())
    set_b = set(splits_b.keys())
    shared = set_a & set_b
    diff = set_a ^ set_b
    max_rf = 2 * (n_tips - 3)
    return {
        "rf": len(diff),
        "max_rf": max_rf,
        "nrf": round(len(diff) / max_rf, 3),
        "shared": len(shared),
        "unique_a": len(set_a - set_b),
        "unique_b": len(set_b - set_a),
    }

# ------------------------------------------------------------------------------
# 4. CURATED DEEP CLADES & TAXON MONOPHYLY
# ------------------------------------------------------------------------------
CURATED_CLADES = [
    ("Bahamas", {"biminiensis", "distichoides", "distichus", "dapsilis", "ocior"}),
    ("Bahamas_minus_ocior", {"biminiensis", "distichoides", "distichus", "dapsilis"}),
    ("(biminiensis, distichoides)", {"biminiensis", "distichoides"}),
    ("((bim, dist), distichus)", {"biminiensis", "distichoides", "distichus"}),
    ("Bahamas + dom_E + ocior", {"biminiensis", "distichoides", "distichus", "dapsilis", "ocior", "dominicensis_E"}),
    ("dominicensis_E", {"dominicensis_E"}),
    ("Western clade (excl. aurifer)", {"dominicensis_NW", "dominicensis_SW", "dominicensis_W", "favillarum", "suppar", "vinosus"}),
    ("(suppar, vinosus)", {"suppar", "vinosus"}),
    ("East complex", {"east_ignigularis", "properus", "sejunctus"}),
    ("(properus, sejunctus)", {"properus", "sejunctus"}),
    ("(east_ignigularis, properus)", {"east_ignigularis", "properus"}),
    ("(west_ignigularis, ravitergum)", {"west_ignigularis", "ravitergum"}),
    ("(dominicensis_NW, dominicensis_W)", {"dominicensis_NW", "dominicensis_W"}),
]


def clade_stats(info, ids):
    cands = [nid for nid, d in info.items() if d["leafset"].issuperset(ids)]
    mn = min(cands, key=lambda nid: len(info[nid]["leafset"]))
    mono = info[mn]["leafset"] == ids
    core_cands = [nid for nid, d in info.items() if d["leafset"].issubset(ids) and len(d["leafset"]) > 1]
    if core_cands:
        core = max(core_cands, key=lambda nid: len(info[nid]["leafset"]))
        core_size = len(info[core]["leafset"])
        core_sup = parse_support(info[core]["label"])
    else:
        core_size = 1
        core_sup = None
    return mono, core_size, len(ids) - core_size, core_sup


def monophyly_info(info, per_tax):
    out = {}
    for tax, ids in per_tax.items():
        cands = [nid for nid, d in info.items() if d["leafset"].issuperset(ids)]
        mn = min(cands, key=lambda nid: len(info[nid]["leafset"]))
        mono = info[mn]["leafset"] == ids
        sup = parse_support(info[mn]["label"])
        core_cands = [nid for nid, d in info.items() if d["leafset"].issubset(ids) and len(d["leafset"]) > 1]
        if core_cands:
            core = max(core_cands, key=lambda nid: len(info[nid]["leafset"]))
            core_size = len(info[core]["leafset"])
            core_sup = parse_support(info[core]["label"])
        else:
            core_size = 1
            core_sup = None
        misplaced = (ids - (info[core]["leafset"] if core_cands else frozenset())) if len(ids) > 1 else frozenset()
        out[tax] = {
            "n": len(ids),
            "monophyletic": mono,
            "clade_support": sup,
            "core_size": core_size,
            "core_support": core_sup,
            "misplaced": sorted(misplaced),
        }
    return out


def destination_of_tip(tip, tax, info):
    cands = [nid for nid, d in info.items() if tip in d["leafset"] and not d["taxa"].issubset({tax})]
    mn = min(cands, key=lambda nid: len(info[nid]["leafset"]))
    d = info[mn]
    return {
        "n_tips": len(d["leafset"]),
        "taxa": sorted(d["taxa"]),
        "support": parse_support(d["label"]),
        "leafset": sorted(d["leafset"]),
    }

# ------------------------------------------------------------------------------
# 5. COLLAPSED BACKBONE RF
# ------------------------------------------------------------------------------
def maximal_blocks(t, tip_to_tax):
    info = annotate(t, tip_to_tax)
    blocks = {}
    for nid, d in info.items():
        if len(d["leafset"]) < 2 or len(d["taxa"]) != 1:
            continue
        p = d["parent"]
        if p is not None and len(info[p]["taxa"]) == 1:
            continue
        blocks[nid] = next(iter(d["taxa"]))
    return blocks, info


def contract_tree(t, tip_to_tax):
    blocks, info = maximal_blocks(t, tip_to_tax)
    if not blocks:
        return None, None
    kids = {}
    for nid, d in info.items():
        kids.setdefault(d["parent"], []).append(nid)
    block_of_node = {}
    names = {}
    seen = Counter()
    for root_nid, tax in sorted(blocks.items()):
        seen[tax] += 1
        name = f"{tax}#{seen[tax]}"
        names[root_nid] = name
        stack = [root_nid]
        while stack:
            u = stack.pop()
            block_of_node[u] = name
            stack.extend(kids.get(u, []))
    named = {root_nid: names[root_nid] for root_nid in sorted(blocks, key=lambda n: (blocks[n], n))}

    def leafset_blocks(nid):
        if nid in block_of_node:
            return frozenset([block_of_node[nid]])
        out = frozenset()
        for c in kids.get(nid, []):
            out |= leafset_blocks(c)
        return out

    res = {nid: leafset_blocks(nid) for nid in info}
    return named, res


def split_sets(sets_by_node, common):
    out = set()
    for s in sets_by_node.values():
        inter = s & common
        if len(inter) == 0 or len(inter) == len(common):
            continue
        comp = common - inter
        small = inter if len(inter) <= len(comp) else comp
        if len(small) >= 2:
            out.add(frozenset(small))
    return out


def compute_backbone_rf(t1, t2, tip_to_tax):
    named1, res1 = contract_tree(t1, tip_to_tax)
    named2, res2 = contract_tree(t2, tip_to_tax)
    common = set(named1.values()) & set(named2.values())
    s1 = split_sets(res1, common)
    s2 = split_sets(res2, common)
    diff = s1 ^ s2
    shared = s1 & s2
    max_rf = 2 * (len(common) - 3) if len(common) >= 3 else 0
    return {
        "n_blocks_1": len(named1),
        "n_blocks_2": len(named2),
        "shared_blocks": len(common),
        "rf": len(diff),
        "max_rf": max_rf,
        "nrf": round(len(diff) / max_rf, 3) if max_rf > 0 else 0.0,
        "shared_splits": len(shared),
        "unique_1": len(s1 - s2),
        "unique_2": len(s2 - s1),
    }

# ------------------------------------------------------------------------------
# 6. FORMATTING HELPERS
# ------------------------------------------------------------------------------
def fmt_sup(sup):
    if not sup:
        return "."
    if sup[1] is None:
        return f"{sup[0]:g}"
    return f"{sup[0]:g}/{sup[1]:g}"

# ------------------------------------------------------------------------------
# 7. MAIN EXECUTION PIPELINE
# ------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Distichus Chromosomal Discordance Analysis Pipeline")
    parser.add_argument("--skip-figures", action="store_true", help="Skip rendering R tanglegram figures")
    parser.add_argument("--outdir", default=None, help="Directory to save generated tables and reports")
    args = parser.parse_args()

    # Determine script location and search paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.abspath(os.path.join(script_dir, ".."))
    outdir = os.path.abspath(args.outdir) if args.outdir else os.path.join(script_dir, "output")
    os.makedirs(outdir, exist_ok=True)

    print("=" * 80)
    print("DISTICHUS CHROMOSOMAL DISCORDANCE ANALYSIS PIPELINE")
    print("=" * 80)
    print(f"Working Directory: {os.getcwd()}")
    print(f"Output Directory:  {outdir}\n")

    # Locate input files
    tree_candidates = {
        "autosomes": [
            os.path.join(project_dir, "trees", "distichus_male_90complete_autosomes.treefile"),
            os.path.join(script_dir, "distichus_male_90complete_autosomes.treefile"),
            os.path.join(project_dir, "trees", "analysis", "distichus_male_90complete_autosomes.treefile"),
        ],
        "Y": [
            os.path.join(project_dir, "trees", "y_distichus_50complete.treefile"),
            os.path.join(script_dir, "y_distichus_50complete.treefile"),
            os.path.join(project_dir, "trees", "analysis", "y_distichus_50complete.treefile"),
        ],
        "X1": [
            os.path.join(project_dir, "trees", "x1_distichus_50complete.treefile"),
            os.path.join(script_dir, "x1_distichus_50complete.treefile"),
            os.path.join(project_dir, "trees", "analysis", "x1_distichus_50complete.treefile"),
        ],
        "X2": [
            os.path.join(project_dir, "trees", "x2_distichus_50complete.treefile"),
            os.path.join(script_dir, "x2_distichus_50complete.treefile"),
            os.path.join(project_dir, "trees", "analysis", "x2_distichus_50complete.treefile"),
        ],
    }

    tree_paths = {}
    for k, cands in tree_candidates.items():
        found = find_file(cands)
        if not found:
            sys.exit(f"ERROR: Could not locate tree file for {k}. Please check your tree directory.")
        tree_paths[k] = found
        print(f"  {k:10} tree: {found}")

    # Metadata files
    tt_file = find_file([
        os.path.join(project_dir, "trees", "analysis", "tips_taxon.csv"),
        os.path.join(script_dir, "tips_taxon.csv"),
    ])
    if not tt_file:
        sys.exit("ERROR: Could not locate tips_taxon.csv.")
    print(f"  Metadata:         {tt_file}")

    loc_file = find_file([
        os.path.join(project_dir, "trees", "analysis", "tip_locality.csv"),
        os.path.join(script_dir, "tip_locality.csv"),
    ])
    print(f"  Locality data:    {loc_file or '(none)'}\n")

    # Load metadata
    tip_to_tax = {}
    with open(tt_file) as fh:
        for row in csv.DictReader(fh):
            tip_to_tax[row["tip"]] = row["taxon"]

    locality = {}
    if loc_file:
        with open(loc_file) as fh:
            for row in csv.DictReader(fh):
                locality[row["tip"]] = row["locality"]

    # Parse raw trees
    raw_trees = {}
    for k, p in tree_paths.items():
        with open(p) as fh:
            raw_trees[k] = parse_newick(fh.read().strip())

    all_tips = frozenset(norm_leaf(l) for l in get_leaves(raw_trees["autosomes"]))
    for k, tr in raw_trees.items():
        assert frozenset(norm_leaf(l) for l in get_leaves(tr)) == all_tips, f"Tip set mismatch in {k}"
    print(f"All trees successfully validated across {len(all_tips)} male specimens (18 taxa).\n")

    # Root configurations
    seeds_auto_x = ("5594 5593 5586 5683 7459 7458 4329 4250 4242 4257 4041 4276 4277 "
                    "5100 5099 5101 590 5102 5094 606 "
                    "9840 9839 9837 9833 9832 4549 4548 4428 4384 4382 4030 7697 3705 "
                    "3703").split()

    seeds_y = ("6769 6766 6768 6617 6789 6787 6453 748 6685_rep "
               "5594 5593 5586 5683 7459 7458 4329 4250 4242 4257 4041 4276 4277 "
               "5100 5099 5101 590 5102 5094 606 "
               "9840 9839 9837 9833 9832 4549 4548 4428 4384 4382 4030 7697 3705 "
               "3703").split()

    basal_seeds = {
        "autosomes": seeds_auto_x,
        "Y": seeds_y,
        "X1": seeds_auto_x,
        "X2": seeds_auto_x,
    }

    # Root each tree and extract annotated node info
    rooted_trees = {}
    annotated_infos = {}
    splits_by_tree = {}
    tree_stats = {}

    for k, tr in raw_trees.items():
        nid, _info, nodes = smallest_clade_node(tr, basal_seeds[k], tip_to_tax)
        r_tree = reroot_edge(tr, nodes[nid])
        rooted_trees[k] = r_tree
        info = annotate(r_tree, tip_to_tax)
        annotated_infos[k] = info
        sp = splits_with_support(info, all_tips, tip_to_tax)
        splits_by_tree[k] = sp
        tree_stats[k] = support_stats(sp)

    # --------------------------------------------------------------------------
    # TABLE 2: SUPPORT STATISTICS
    # --------------------------------------------------------------------------
    print("=" * 80)
    print("TABLE 2. TREE SUPPORT STATISTICS")
    print("=" * 80)
    tbl2_rows = [
        ("Tips", {k: len(all_tips) for k in tree_stats}),
        ("Taxa", {k: 18 for k in tree_stats}),
        ("Total splits", {k: tree_stats[k]["n_splits"] for k in tree_stats}),
        ("Mean SH-aLRT", {k: f"{tree_stats[k]['shalrt_mean']:.1f}" for k in tree_stats}),
        ("Median SH-aLRT", {k: f"{tree_stats[k]['shalrt_median']:.1f}" for k in tree_stats}),
        ("Mean UFBoot", {k: f"{tree_stats[k]['ufboot_mean']:.1f}" for k in tree_stats}),
        ("Median UFBoot", {k: f"{tree_stats[k]['ufboot_median']:.1f}" for k in tree_stats}),
        ("% SH-aLRT >= 95", {k: f"{tree_stats[k]['pct_shalrt_ge95']:.1f}%" for k in tree_stats}),
        ("% UFBoot >= 95", {k: f"{tree_stats[k]['pct_ufboot_ge95']:.1f}%" for k in tree_stats}),
        ("% SH-aLRT = 100", {k: f"{tree_stats[k]['pct_both100']:.1f}%" for k in tree_stats}),
    ]

    header_t2 = f"{'Metric':<22} | {'Autosomes':<12} | {'Y Chromosome':<12} | {'X1 (Ancestral)':<14} | {'X2 (Neo-X)':<12}"
    print(header_t2)
    print("-" * len(header_t2))
    for m, vals in tbl2_rows:
        print(f"{m:<22} | {vals['autosomes']:<12} | {vals['Y']:<12} | {vals['X1']:<14} | {vals['X2']:<12}")
    print()

    # Save Table 2 CSV
    with open(os.path.join(outdir, "table2_support_statistics.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Metric", "Autosomes", "Y_Chromosome", "X1_Ancestral", "X2_NeoX"])
        for m, vals in tbl2_rows:
            w.writerow([m, vals["autosomes"], vals["Y"], vals["X1"], vals["X2"]])

    # --------------------------------------------------------------------------
    # TABLE 3: ROBINSON-FOULDS DISTANCES
    # --------------------------------------------------------------------------
    print("=" * 80)
    print("TABLE 3. ROBINSON-FOULDS DISTANCE COMPARISONS")
    print("=" * 80)
    pairs = [
        ("Autosomes", "Y", splits_by_tree["autosomes"], splits_by_tree["Y"]),
        ("X1", "X2", splits_by_tree["X1"], splits_by_tree["X2"]),
        ("Autosomes", "X1", splits_by_tree["autosomes"], splits_by_tree["X1"]),
        ("Autosomes", "X2", splits_by_tree["autosomes"], splits_by_tree["X2"]),
        ("Y", "X1", splits_by_tree["Y"], splits_by_tree["X1"]),
        ("Y", "X2", splits_by_tree["Y"], splits_by_tree["X2"]),
    ]
    rf_results = {}
    header_t3 = f"{'Comparison':<20} | {'Raw RF':<8} | {'Max RF':<8} | {'nRF':<8} | {'Shared Splits':<14} | {'Unique (Tree 1)':<15} | {'Unique (Tree 2)':<15}"
    print(header_t3)
    print("-" * len(header_t3))
    with open(os.path.join(outdir, "table3_robinson_foulds.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Comparison", "Raw_RF", "Max_RF", "nRF", "Shared_Splits", "Unique_Tree1", "Unique_Tree2"])
        for n1, n2, s1, s2 in pairs:
            res = compute_rf(s1, s2)
            rf_results[f"{n1}_vs_{n2}"] = res
            cmp_name = f"{n1} vs. {n2}"
            print(f"{cmp_name:<20} | {res['rf']:<8} | {res['max_rf']:<8} | {res['nrf']:<8.3f} | {res['shared']:<14} | {res['unique_a']:<15} | {res['unique_b']:<15}")
            w.writerow([cmp_name, res["rf"], res["max_rf"], res["nrf"], res["shared"], res["unique_a"], res["unique_b"]])
    print()

    # --------------------------------------------------------------------------
    # TABLE 4: CURATED DEEP CLADE MONOPHYLY AND SUPPORT
    # --------------------------------------------------------------------------
    print("=" * 80)
    print("TABLE 4. DEEP CLADE MONOPHYLY AND SUPPORT ACROSS ALL GENOMIC COMPARTMENTS")
    print("=" * 80)
    pt = defaultdict(set)
    for tip, tax in tip_to_tax.items():
        pt[tax].add(tip)

    deep_clade_results = []
    header_t4 = f"{'Clade':<28} | {'n':<3} | {'Auto mono (core, sup)':<22} | {'Y mono (core, sup)':<22} | {'X1 mono (core, sup)':<22} | {'X2 mono (core, sup)':<22}"
    print(header_t4)
    print("-" * len(header_t4))
    with open(os.path.join(outdir, "table4_deep_clades.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Clade", "n_tips",
                    "Auto_mono", "Auto_core", "Auto_support",
                    "Y_mono", "Y_core", "Y_support",
                    "X1_mono", "X1_core", "X1_support",
                    "X2_mono", "X2_core", "X2_support"])
        for name, taxa in CURATED_CLADES:
            ids = set()
            for tx in taxa:
                ids |= pt[tx]
            row_dict = {"name": name, "n": len(ids)}
            row_print = [f"{name:<28}", f"{len(ids):<3}"]
            csv_row = [name, len(ids)]
            for k in ["autosomes", "Y", "X1", "X2"]:
                m, c, miss, sup = clade_stats(annotated_infos[k], ids)
                sup_s = fmt_sup(sup)
                row_dict[k] = {"mono": m, "core": c, "sup": sup_s}
                row_print.append(f"{'Yes' if m else 'No':<3} ({c:<2}, {sup_s:<9})")
                csv_row.extend(["Yes" if m else "No", c, sup_s])
            print(" | ".join(row_print))
            w.writerow(csv_row)
            deep_clade_results.append(row_dict)
    print()

    # --------------------------------------------------------------------------
    # TABLE 5: TAXON-LEVEL MONOPHYLY & DISCORDANCE
    # --------------------------------------------------------------------------
    print("=" * 80)
    print("TABLE 5. TAXON-LEVEL MONOPHYLY & STRAYS")
    print("=" * 80)
    mono_all = {k: monophyly_info(annotated_infos[k], pt) for k in annotated_infos}

    # Taxon table for Autosomes vs Y and X1 vs X2
    with open(os.path.join(outdir, "table5_taxon_discordance_complete.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "Taxon", "n",
            "Auto_core", "Auto_strays", "Y_core", "Y_strays", "Auto_Y_Class",
            "X1_core", "X1_strays", "X2_core", "X2_strays", "X1_X2_Class"
        ])
        header_t5 = f"{'Taxon':<18} | {'n':<3} | {'Auto Y-Class':<18} | {'(Auto/Y strays)':<16} | {'X1 X2-Class':<22} | {'(X1/X2 strays)':<16}"
        print(header_t5)
        print("-" * len(header_t5))
        for tax in sorted(pt):
            n_t = len(pt[tax])
            # Auto vs Y
            ma = mono_all["autosomes"][tax]
            my = mono_all["Y"][tax]
            sa = len(ma["misplaced"])
            sy = len(my["misplaced"])
            cls_ay = "Concordant" if sa == 0 and sy == 0 else ("Discordant (Auto)" if sy == 0 else ("Discordant (Y)" if sa == 0 else "Discordant (Both)"))

            # X1 vs X2
            m1 = mono_all["X1"][tax]
            m2 = mono_all["X2"][tax]
            s1 = len(m1["misplaced"])
            s2 = len(m2["misplaced"])
            cls_12 = "Concordant" if s1 == 0 and s2 == 0 else ("Discordant (X1)" if s2 == 0 else ("Discordant (X2)" if s1 == 0 else "Discordant (Both)"))

            print(f"{tax:<18} | {n_t:<3} | {cls_ay:<18} | {f'{sa} / {sy}':<16} | {cls_12:<22} | {f'{s1} / {s2}':<16}")
            w.writerow([
                tax, n_t,
                ma["core_size"], sa, my["core_size"], sy, cls_ay,
                m1["core_size"], s1, m2["core_size"], s2, cls_12
            ])
    print()

    # --------------------------------------------------------------------------
    # TABLE 6: INDIVIDUAL-LEVEL DISCORDANCE & STRAYS
    # --------------------------------------------------------------------------
    print("=" * 80)
    print("TABLE 6. INDIVIDUAL-LEVEL DISCORDANCE BREAKDOWN")
    print("=" * 80)
    strays_by_tree = {}
    for k in mono_all:
        strays_by_tree[k] = {tip: tax for tax, m in mono_all[k].items() for tip in m["misplaced"]}

    def breakdown_pair(k1, k2):
        s1 = set(strays_by_tree[k1].keys())
        s2 = set(strays_by_tree[k2].keys())
        both = s1 & s2
        only1 = s1 - s2
        only2 = s2 - s1
        total_strays = s1 | s2
        concordant = len(all_tips) - len(total_strays)
        return {
            "concordant": (concordant, round(100 * concordant / len(all_tips), 1)),
            f"only_{k1}": (len(only1), round(100 * len(only1) / len(all_tips), 1)),
            f"only_{k2}": (len(only2), round(100 * len(only2) / len(all_tips), 1)),
            "both": (len(both), round(100 * len(both) / len(all_tips), 1)),
            "total_strays": (len(total_strays), round(100 * len(total_strays) / len(all_tips), 1)),
        }

    bd_ay = breakdown_pair("autosomes", "Y")
    bd_12 = breakdown_pair("X1", "X2")

    ay_conc = f"{bd_ay['concordant'][0]} ({bd_ay['concordant'][1]}%)"
    x_conc  = f"{bd_12['concordant'][0]} ({bd_12['concordant'][1]}%)"
    ay_t1   = f"{bd_ay['only_autosomes'][0]} ({bd_ay['only_autosomes'][1]}%)"
    x_t1    = f"{bd_12['only_X1'][0]} ({bd_12['only_X1'][1]}%)"
    ay_t2   = f"{bd_ay['only_Y'][0]} ({bd_ay['only_Y'][1]}%)"
    x_t2    = f"{bd_12['only_X2'][0]} ({bd_12['only_X2'][1]}%)"
    ay_both = f"{bd_ay['both'][0]} ({bd_ay['both'][1]}%)"
    x_both  = f"{bd_12['both'][0]} ({bd_12['both'][1]}%)"
    ay_tot  = f"{bd_ay['total_strays'][0]} ({bd_ay['total_strays'][1]}%)"
    x_tot   = f"{bd_12['total_strays'][0]} ({bd_12['total_strays'][1]}%)"

    print(f"{'Category':<28} | {'Autosomes vs. Y (n, %)':<25} | {'X1 vs. X2 (n, %)':<25}")
    print("-" * 82)
    print(f"{'Concordant in both':<28} | {ay_conc:<25} | {x_conc:<25}")
    print(f"{'Stray in Tree 1 only':<28} | {ay_t1:<25} | {x_t1:<25}")
    print(f"{'Stray in Tree 2 only':<28} | {ay_t2:<25} | {x_t2:<25}")
    print(f"{'Stray in both genomes':<28} | {ay_both:<25} | {x_both:<25}")
    print(f"{'Total strays (>= 1 genome)':<28} | {ay_tot:<25} | {x_tot:<25}")
    print()

    # Detailed X1 vs X2 stray tip catalog
    with open(os.path.join(outdir, "individual_strays_x1_x2.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Tip", "Taxon", "Locality", "X1_Stray", "X1_Dest_Taxa", "X1_Dest_n", "X1_Dest_Sup",
                    "X2_Stray", "X2_Dest_Taxa", "X2_Dest_n", "X2_Dest_Sup", "Classification"])
        for tip in sorted(all_tips):
            in_x1 = tip in strays_by_tree["X1"]
            in_x2 = tip in strays_by_tree["X2"]
            tax = tip_to_tax[tip]
            loc = locality.get(tip, "?")
            cls = "concordant" if (not in_x1 and not in_x2) else ("both" if in_x1 and in_x2 else ("X1_only" if in_x1 else "X2_only"))

            d1 = destination_of_tip(tip, tax, annotated_infos["X1"]) if in_x1 else None
            d2 = destination_of_tip(tip, tax, annotated_infos["X2"]) if in_x2 else None

            w.writerow([
                tip, tax, loc,
                "1" if in_x1 else "0", "+".join(d1["taxa"]) if d1 else ".", d1["n_tips"] if d1 else "0", fmt_sup(d1["support"]) if d1 else ".",
                "1" if in_x2 else "0", "+".join(d2["taxa"]) if d2 else ".", d2["n_tips"] if d2 else "0", fmt_sup(d2["support"]) if d2 else ".",
                cls
            ])

    # --------------------------------------------------------------------------
    # 7. COLLAPSED BACKBONE RF ANALYSIS
    # --------------------------------------------------------------------------
    print("=" * 80)
    print("COLLAPSED BACKBONE ROBINSON-FOULDS ANALYSIS")
    print("=" * 80)
    bb_ay = compute_backbone_rf(rooted_trees["autosomes"], rooted_trees["Y"], tip_to_tax)
    bb_12 = compute_backbone_rf(rooted_trees["X1"], rooted_trees["X2"], tip_to_tax)

    print(f"Autosomes vs. Y:")
    print(f"  Maximal blocks:  25 Autosomal vs. 21 Y-chromosome (Comparable shared blocks: {bb_ay['shared_blocks']})")
    print(f"  Backbone RF:     {bb_ay['rf']} / {bb_ay['max_rf']} (nRF = {bb_ay['nrf']:.3f})")
    print(f"  Backbone Splits: {bb_ay['shared_splits']} shared, {bb_ay['unique_1']} unique to Autosomes, {bb_ay['unique_2']} unique to Y\n")

    print(f"Ancestral X (X1) vs. Neo-X (X2):")
    print(f"  Maximal blocks:  {bb_12['n_blocks_1']} X1 vs. {bb_12['n_blocks_2']} X2 (Comparable shared blocks: {bb_12['shared_blocks']})")
    print(f"  Backbone RF:     {bb_12['rf']} / {bb_12['max_rf']} (nRF = {bb_12['nrf']:.3f})")
    print(f"  Backbone Splits: {bb_12['shared_splits']} shared, {bb_12['unique_1']} unique to X1, {bb_12['unique_2']} unique to X2\n")

    # --------------------------------------------------------------------------
    # 8. FIGURE GENERATION (VIA RSCRIPT)
    # --------------------------------------------------------------------------
    if not args.skip_figures:
        print("=" * 80)
        print("GENERATING PUBLICATION TANGLEGRAM FIGURES")
        print("=" * 80)

        # 0. Figure 1: Multi-chromosomal discordance at three resolutions
        py_fig1 = find_file([
            os.path.join(script_dir, "plot_discordance_levels.py"),
            os.path.join(project_dir, "trees", "analysis", "plot_discordance_levels.py"),
        ])
        if py_fig1:
            print(f"Running Figure 1 renderer ({py_fig1})...")
            try:
                subprocess.run([sys.executable, py_fig1], cwd=script_dir, check=True)
                print("  -> Generated discordance_levels.pdf and discordance_levels.png successfully.")
            except Exception as e:
                print(f"  WARNING: Could not run Figure 1 renderer: {e}")

        # 1. Autosomes vs Y
        r_script_ay = find_file([
            os.path.join(script_dir, "distichus_autosomes_Y_discordance_tanglegram.R"),
            os.path.join(project_dir, "trees", "analysis", "plot_discordance_tanglegram_fig2_style.R"),
        ])
        if r_script_ay:
            print(f"Running Autosomes vs. Y tanglegram renderer ({r_script_ay})...")
            try:
                subprocess.run(["Rscript", r_script_ay], cwd=script_dir, check=True)
                print("  -> Generated discordance_tanglegram_Y_auto.pdf successfully.")
            except Exception as e:
                print(f"  WARNING: Could not run Rscript for Autosomes vs Y: {e}")

        # 2. X1 vs X2
        r_script_12 = find_file([
            os.path.join(script_dir, "distichus_x1_x2_tanglegram.R"),
            os.path.join(project_dir, "trees", "analysis", "distichus_x1_x2_tanglegram.R"),
        ])
        if r_script_12:
            print(f"Running X1 vs. X2 tanglegram renderer ({r_script_12})...")
            try:
                subprocess.run(["Rscript", r_script_12], cwd=script_dir, check=True)
                print("  -> Generated discordance_tanglegram_x1_x2.pdf successfully.")
            except Exception as e:
                print(f"  WARNING: Could not run Rscript for X1 vs X2: {e}")
        print()

    # Save summary JSON
    summary_payload = {
        "support_statistics": tree_stats,
        "robinson_foulds": rf_results,
        "curated_deep_clades": deep_clade_results,
        "backbone_rf": {
            "autosomes_vs_y": bb_ay,
            "x1_vs_x2": bb_12,
        },
        "individual_breakdown": {
            "autosomes_vs_y": bd_ay,
            "x1_vs_x2": bd_12,
        },
    }
    with open(os.path.join(outdir, "discordance_complete_summary.json"), "w") as fh:
        json.dump(summary_payload, fh, indent=2)

    print("=" * 80)
    print("ANALYSIS COMPLETE!")
    print(f"Summary data and publication tables written to: {outdir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
