#!/usr/bin/env python3
"""Emit canonical JSON results from the Python taxutils implementation."""

import json
import math
import sys

import numpy as np
import pandas as pd

from taxutils import taxutils


def native(value):
    """Convert pandas/NumPy values into stable JSON-compatible values."""
    if isinstance(value, np.generic):
        value = value.item()
    if value is pd.NA or (isinstance(value, float) and math.isnan(value)):
        return None
    if isinstance(value, dict):
        return {str(key): native(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray, pd.Series)):
        return [native(item) for item in list(value)]
    return value


def profile(series):
    return {str(key): native(value) for key, value in series.to_dict().items()}


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: python_reference.py FIXTURE_DIR")

    tu = taxutils(low_memory=True)
    accession_inputs = [
        ">NC_000001.1 bacterial alpha",
        ">kraken:taxid|22|ABCD01000001.1 viral example",
        "prefix NC_000003.2 suffix",
        "nc_000002.1",
        "sequence_without_accession",
        "FOONC_045512.2",
    ]
    node_columns = [
        "taxon",
        "parent",
        "rank",
        "rank_code",
        "rank_base",
        "rank_idx",
        "new_rank",
    ]

    tu.load_a2t(
        [
            ">NC_000001.1 description",
            "ABCD01000001.1",
            "missing",
        ]
    )
    a2t_replace = sorted((str(key), native(value)) for key, value in tu.a2t.items())
    tu.load_a2t(["NC_000002.1", "NC_000001.1"], extend=True)
    a2t_extend = sorted((str(key), native(value)) for key, value in tu.a2t.items())
    tu.load_a2t(["WXYZ01000001.1"], wgs=True)
    a2t_wgs = sorted((str(key), native(value)) for key, value in tu.a2t.items())
    tu.load_a2t(
        ["NC_000001.1", "WXYZ01000001.1"],
        low_memory=False,
        wgs=True,
    )
    a2t_sqlite_wgs = sorted(
        (str(key), native(value)) for key, value in tu.a2t.items()
    )

    tree_with_ancestors = tu.format_tree([14, 15, 22])
    tree_without_ancestors = tu.format_tree(
        [14, 15, 22], include_ancestors=False
    )
    topology_many = tu.topology([13, 15], anchor_rank="F")
    topology_stats = [
        "n_taxa",
        "n_leaves",
        "max_depth",
        "mean_depth",
        "topology_scale",
        "max_children",
        "branching_taxa_fraction",
        "top_child_fraction",
    ]

    result = {
        "accessions_versioned": native(tu.parse_accession(accession_inputs)),
        "accessions_unversioned": native(
            tu.parse_accession(accession_inputs, version=False)
        ),
        "rank_order": tu.get_rank_order(),
        "names": [
            [taxon, tu.names[taxon]]
            for taxon in [0, 1, 2, 13, 22, 694009, 2697049]
        ],
        "nodes": native(tu.nodes[node_columns].to_dict("records")),
        "parents": [
            [taxon, native(tu.parent.get(taxon))]
            for taxon in [1, 2, 3, 10, 11, 12, 13, 14, 15, 10239, 20, 21, 22]
        ],
        "target_taxa": native(tu.target_taxa),
        "branch": native(tu.get_branch(14)),
        "subtree": native(tu.get_subtree(12)),
        "ancestors": native(tu.get_ancestor([14, 15, 22], "F")),
        "leaf_scalar": native(tu.is_leaf(14)),
        "leaf_many": native(tu.is_leaf([12, 14, 15, 999999])),
        "child_scalar": native(tu.is_child(14, 13)),
        "child_many": native(tu.is_child([14, 15, 22], [13, 12, 21])),
        "descendent_scalar": native(tu.is_descendent(14, 10)),
        "descendent_many": native(
            tu.is_descendent([14, 13, 10, 22], [10, 13, 10, 20])
        ),
        "lca": native(tu.get_lca(14, 15)),
        "distance": native(tu.get_distance(14, 15)),
        "sort_taxa": native(tu.sort_taxa([22, 15, 10, 14, 13, 15])),
        "tree_with_ancestors": [
            [int(taxon), name] for taxon, name in tree_with_ancestors.items()
        ],
        "tree_without_ancestors": [
            [int(taxon), name] for taxon, name in tree_without_ancestors.items()
        ],
        "higher_than_family": native(
            tu.higher_than_rank([1, 2, 10, 11, 12, 13, 999999], "F")
        ),
        "topology": profile(tu.topology(12)),
        "topology_family_anchor": profile(tu.topology(14, anchor_rank="F")),
        "topology_many": [
            {str(key): native(value) for key, value in row.items()}
            for row in topology_many.to_dict("records")
        ],
        "topology_scale": native(
            tu.topology([13, 15], anchor_rank="F", stat="topology_scale")
        ),
        "topology_stats": {
            stat: float(tu.topology(12, stat=stat))
            for stat in topology_stats
        },
        "a2t_replace": a2t_replace,
        "a2t_extend": a2t_extend,
        "a2t_wgs": a2t_wgs,
        "a2t_sqlite_wgs": a2t_sqlite_wgs,
        "t2a": sorted(tu.get_t2a([13, 22])),
        "t2a_sqlite_wgs": sorted(
            tu.get_t2a([15, 22], low_memory=False, wgs=True)
        ),
    }
    json.dump(result, sys.stdout, sort_keys=True, separators=(",", ":"))


if __name__ == "__main__":
    main()
