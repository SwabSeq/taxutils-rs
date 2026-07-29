use std::collections::{BTreeSet, HashMap, HashSet};

use anyhow::{Result, bail};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::resources::{TaxutilsOptions, load_taxutils};

pub type TaxonId = i64;
pub(crate) const RANK_CODES: [&str; 10] = ["U", "R", "D", "K", "P", "C", "O", "F", "G", "S"];

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct TaxonNode {
    pub taxon: TaxonId,
    pub parent: Option<TaxonId>,
    /// Original, lower-case NCBI rank.
    pub rank: String,
    /// Corrected rank code, including subranks (`F2`, `S3`, ...).
    pub rank_code: String,
    /// Canonical base rank code.
    pub rank_base: char,
    /// Canonical rank order index.
    pub rank_idx: u8,
    /// Corrected canonical rank name.
    pub new_rank: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TopologyStat {
    NTaxa,
    NLeaves,
    MaxDepth,
    MeanDepth,
    TopologyScale,
    MaxChildren,
    BranchingTaxaFraction,
    TopChildFraction,
}

impl std::str::FromStr for TopologyStat {
    type Err = anyhow::Error;
    fn from_str(value: &str) -> Result<Self> {
        match value {
            "n_taxa" => Ok(Self::NTaxa),
            "n_leaves" => Ok(Self::NLeaves),
            "max_depth" => Ok(Self::MaxDepth),
            "mean_depth" => Ok(Self::MeanDepth),
            "topology_scale" => Ok(Self::TopologyScale),
            "max_children" => Ok(Self::MaxChildren),
            "branching_taxa_fraction" => Ok(Self::BranchingTaxaFraction),
            "top_child_fraction" => Ok(Self::TopChildFraction),
            _ => bail!(
                "stat must be one of: n_taxa, n_leaves, max_depth, mean_depth, \
                 topology_scale, max_children, branching_taxa_fraction, top_child_fraction"
            ),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TopologyProfile {
    pub taxon: TaxonId,
    pub name: String,
    pub rank_code: Option<String>,
    pub anchor_taxon: TaxonId,
    pub anchor_name: String,
    pub anchor_rank_code: Option<String>,
    pub n_taxa: usize,
    pub n_leaves: usize,
    pub max_depth: usize,
    pub mean_depth: f64,
    pub topology_scale: usize,
    pub max_children: usize,
    pub branching_taxa_fraction: f64,
    pub top_child_fraction: f64,
}

#[derive(Debug)]
pub struct TaxonomicUtils {
    pub names: HashMap<TaxonId, String>,
    pub nodes: Vec<TaxonNode>,
    pub target_taxa: Vec<TaxonId>,
    pub a2t: HashMap<String, TaxonId>,
    pub parent: HashMap<TaxonId, Option<TaxonId>>,
    pub(crate) low_memory: bool,
    pub(crate) wgs: bool,
    pub(crate) save_folder: std::path::PathBuf,
    children: HashMap<TaxonId, Vec<TaxonId>>,
    node_index: HashMap<TaxonId, usize>,
    depth: HashMap<TaxonId, usize>,
    descendant_index: HashMap<TaxonId, (usize, usize)>,
}

impl TaxonomicUtils {
    pub fn new(options: TaxutilsOptions) -> Result<Self> {
        load_taxutils(options)
    }

    pub(crate) fn from_parts(
        names: HashMap<TaxonId, String>,
        nodes: Vec<TaxonNode>,
        target_taxa: Vec<TaxonId>,
        a2t: HashMap<String, TaxonId>,
        low_memory: bool,
        wgs: bool,
        save_folder: std::path::PathBuf,
    ) -> Self {
        let mut parent = nodes
            .iter()
            .map(|n| (n.taxon, n.parent))
            .collect::<HashMap<_, _>>();
        parent.insert(1, None);
        let mut children: HashMap<TaxonId, Vec<TaxonId>> = HashMap::new();
        for node in &nodes {
            if node.taxon != 1
                && let Some(parent) = parent.get(&node.taxon).copied().flatten()
            {
                children.entry(parent).or_default().push(node.taxon);
            }
        }
        let node_index = nodes
            .iter()
            .enumerate()
            .map(|(i, node)| (node.taxon, i))
            .collect();
        let (depth, descendant_index) = build_indexes(&parent, &children);
        Self {
            names,
            nodes,
            target_taxa,
            a2t,
            parent,
            low_memory,
            wgs,
            save_folder,
            children,
            node_index,
            depth,
            descendant_index,
        }
    }

    pub fn parse_accession(&self, text: &str, version: bool) -> String {
        crate::accession::parse_accession(text, version)
    }

    pub fn parse_accessions<S: AsRef<str> + Sync>(
        &self,
        values: &[S],
        version: bool,
    ) -> Vec<String> {
        crate::accession::parse_accessions(values, version)
    }

    pub fn get_rank_order(&self) -> Vec<&'static str> {
        RANK_CODES.to_vec()
    }

    pub(crate) fn node(&self, taxon: TaxonId) -> Option<&TaxonNode> {
        self.node_index.get(&taxon).map(|index| &self.nodes[*index])
    }

    pub fn get_branch(&self, taxon: TaxonId) -> Vec<TaxonId> {
        let mut branch = Vec::new();
        let mut current = Some(taxon);
        let mut seen = HashSet::new();
        while let Some(node) = current {
            if !seen.insert(node) {
                break;
            }
            branch.push(node);
            current = self.parent.get(&node).copied().flatten();
        }
        branch.reverse();
        branch
    }

    pub fn get_subtree(&self, taxon: TaxonId) -> Vec<TaxonId> {
        let mut result = Vec::new();
        let mut stack = vec![taxon];
        while let Some(node) = stack.pop() {
            result.push(node);
            if let Some(children) = self.children.get(&node) {
                stack.extend(children.iter().rev().copied());
            }
        }
        result
    }

    pub fn get_ancestor(&self, taxon: TaxonId, anchor_rank: &str) -> Result<TaxonId> {
        let code = rank_to_code(anchor_rank)?;
        for ancestor in self.get_branch(taxon).into_iter().rev() {
            if self.node(ancestor).is_some_and(|n| n.rank_base == code) {
                return Ok(ancestor);
            }
        }
        Ok(taxon)
    }

    pub fn get_ancestors(&self, taxa: &[TaxonId], anchor_rank: &str) -> Result<Vec<TaxonId>> {
        taxa.par_iter()
            .map(|taxon| self.get_ancestor(*taxon, anchor_rank))
            .collect()
    }

    /// Compute topology profiles concurrently while preserving input order.
    pub fn topologies(
        &self,
        taxa: &[TaxonId],
        anchor_rank: Option<&str>,
    ) -> Result<Vec<TopologyProfile>> {
        taxa.par_iter()
            .map(|taxon| self.topology(*taxon, anchor_rank))
            .collect()
    }

    /// Batch form of [`Self::is_leaf`], evaluated in parallel in input order.
    pub fn are_leaves(&self, taxa: &[TaxonId]) -> Vec<bool> {
        taxa.par_iter().map(|taxon| self.is_leaf(*taxon)).collect()
    }

    pub fn is_leaf(&self, taxon: TaxonId) -> bool {
        !self.children.contains_key(&taxon)
    }

    pub fn is_child(&self, taxon_a: TaxonId, taxon_b: TaxonId) -> bool {
        self.parent.get(&taxon_a).copied().flatten() == Some(taxon_b)
    }

    /// Pairwise batch form of [`Self::is_child`], evaluated in parallel.
    pub fn are_children(&self, taxon_a: &[TaxonId], taxon_b: &[TaxonId]) -> Result<Vec<bool>> {
        if taxon_a.len() != taxon_b.len() {
            bail!("taxon_a and taxon_b must have the same length");
        }
        Ok(taxon_a
            .par_iter()
            .zip(taxon_b.par_iter())
            .map(|(a, b)| self.is_child(*a, *b))
            .collect())
    }

    /// Strict descendant check; a taxon is not its own descendant.
    pub fn is_descendent(&self, taxon_a: TaxonId, taxon_b: TaxonId) -> bool {
        if taxon_a == taxon_b {
            return false;
        }
        match (
            self.descendant_index.get(&taxon_a),
            self.descendant_index.get(&taxon_b),
        ) {
            (Some((a, _)), Some((b_start, b_end))) => b_start <= a && a <= b_end,
            _ => false,
        }
    }

    /// Pairwise batch form of [`Self::is_descendent`], evaluated in parallel.
    pub fn are_descendents(&self, taxon_a: &[TaxonId], taxon_b: &[TaxonId]) -> Result<Vec<bool>> {
        if taxon_a.len() != taxon_b.len() {
            bail!("taxon_a and taxon_b must have the same length");
        }
        Ok(taxon_a
            .par_iter()
            .zip(taxon_b.par_iter())
            .map(|(a, b)| self.is_descendent(*a, *b))
            .collect())
    }

    pub fn get_lca(&self, mut a: TaxonId, mut b: TaxonId) -> TaxonId {
        if a == b {
            return a;
        }
        let mut depth_a = self.depth.get(&a).copied().unwrap_or(0);
        let mut depth_b = self.depth.get(&b).copied().unwrap_or(0);
        while depth_a > depth_b {
            let Some(parent) = self.parent.get(&a).copied().flatten() else {
                return 1;
            };
            a = parent;
            depth_a -= 1;
        }
        while depth_b > depth_a {
            let Some(parent) = self.parent.get(&b).copied().flatten() else {
                return 1;
            };
            b = parent;
            depth_b -= 1;
        }
        while a != b {
            let Some(parent_a) = self.parent.get(&a).copied().flatten() else {
                return 1;
            };
            let Some(parent_b) = self.parent.get(&b).copied().flatten() else {
                return 1;
            };
            a = parent_a;
            b = parent_b;
        }
        a
    }

    pub fn get_distance(&self, a: TaxonId, b: TaxonId) -> usize {
        let lca = self.get_lca(a, b);
        self.depth.get(&a).copied().unwrap_or(0) + self.depth.get(&b).copied().unwrap_or(0)
            - 2 * self.depth.get(&lca).copied().unwrap_or(0)
    }

    pub fn higher_than_rank(&self, taxa: &[TaxonId], rank: &str) -> Result<Vec<bool>> {
        let threshold = rank_index(rank_to_code(rank)?);
        Ok(taxa
            .par_iter()
            .map(|taxon| self.node(*taxon).map_or(threshold, |n| n.rank_idx) < threshold)
            .collect())
    }

    pub fn sort_taxa<I>(&self, taxa: I) -> Vec<TaxonId>
    where
        I: IntoIterator<Item = TaxonId>,
    {
        taxonomic_order(
            taxa.into_iter().collect(),
            &self.parent,
            &self
                .nodes
                .iter()
                .map(|n| (n.taxon, n.rank_code.clone()))
                .collect(),
            &self.names,
        )
    }

    pub fn format_tree<I>(
        &self,
        taxa: I,
        include_ancestors: bool,
        root: TaxonId,
        indent: &str,
    ) -> Vec<(TaxonId, String)>
    where
        I: IntoIterator<Item = TaxonId>,
    {
        let taxa: HashSet<_> = taxa.into_iter().collect();
        let mut visible = if include_ancestors {
            HashSet::new()
        } else {
            taxa.clone()
        };
        if include_ancestors {
            for taxon in taxa {
                let mut current = Some(taxon);
                let mut seen = HashSet::new();
                while let Some(node) = current {
                    if !seen.insert(node) {
                        break;
                    }
                    visible.insert(node);
                    if node == root {
                        break;
                    }
                    current = self.parent.get(&node).copied().flatten();
                }
            }
        }
        let order = self.sort_taxa(visible.iter().copied());
        let mut rows = Vec::with_capacity(order.len());
        for taxon in order {
            let depth = self
                .get_branch(taxon)
                .into_iter()
                .filter(|ancestor| visible.contains(ancestor))
                .count()
                .saturating_sub(1);
            rows.push((
                taxon,
                format!(
                    "{}{}",
                    indent.repeat(depth),
                    self.names
                        .get(&taxon)
                        .cloned()
                        .unwrap_or_else(|| taxon.to_string())
                ),
            ));
        }
        rows
    }

    pub fn topology(&self, taxon: TaxonId, anchor_rank: Option<&str>) -> Result<TopologyProfile> {
        let anchor = match anchor_rank {
            Some(rank) => self.get_ancestor(taxon, rank)?,
            None => taxon,
        };
        let subtree = self.get_subtree(anchor);
        let subtree_set: HashSet<_> = subtree.iter().copied().collect();
        let anchor_depth = self.depth.get(&anchor).copied().unwrap_or(0);
        let mut relative_depths = subtree
            .iter()
            .map(|node| {
                self.depth
                    .get(node)
                    .copied()
                    .unwrap_or(0)
                    .saturating_sub(anchor_depth)
            })
            .collect::<Vec<_>>();
        let child_counts = subtree
            .iter()
            .map(|node| {
                self.children.get(node).map_or(0, |values| {
                    values.iter().filter(|v| subtree_set.contains(v)).count()
                })
            })
            .collect::<Vec<_>>();
        let n_taxa = subtree.len();
        let n_leaves = child_counts.iter().filter(|count| **count == 0).count();
        let max_depth = relative_depths.iter().copied().max().unwrap_or(0);
        let mean_depth = if n_taxa == 0 {
            0.0
        } else {
            relative_depths.iter().sum::<usize>() as f64 / n_taxa as f64
        };
        relative_depths.retain(|depth| *depth > 0);
        relative_depths.sort_unstable();
        let topology_scale = relative_depths
            .get((0.95 * relative_depths.len().saturating_sub(1) as f64) as usize)
            .copied()
            .unwrap_or(0)
            .max(1);
        let max_children = child_counts.iter().copied().max().unwrap_or(0);
        let branching_taxa_fraction = if n_taxa == 0 {
            0.0
        } else {
            child_counts.iter().filter(|count| **count > 0).count() as f64 / n_taxa as f64
        };
        let mut subtree_sizes = HashMap::new();
        for node in subtree.iter().rev() {
            let size = 1 + self.children.get(node).map_or(0, |children| {
                children
                    .iter()
                    .filter(|child| subtree_set.contains(child))
                    .map(|child| subtree_sizes.get(child).copied().unwrap_or(0))
                    .sum::<usize>()
            });
            subtree_sizes.insert(*node, size);
        }
        let immediate_sizes = self
            .children
            .get(&anchor)
            .into_iter()
            .flatten()
            .filter(|child| subtree_set.contains(child))
            .map(|child| subtree_sizes.get(child).copied().unwrap_or(0))
            .collect::<Vec<_>>();
        let total_child_size: usize = immediate_sizes.iter().sum();
        let top_child_fraction = if total_child_size == 0 {
            1.0
        } else {
            immediate_sizes.iter().copied().max().unwrap_or(0) as f64 / total_child_size as f64
        };
        Ok(TopologyProfile {
            taxon,
            name: self
                .names
                .get(&taxon)
                .cloned()
                .unwrap_or_else(|| taxon.to_string()),
            rank_code: self.node(taxon).map(|n| n.rank_code.clone()),
            anchor_taxon: anchor,
            anchor_name: self
                .names
                .get(&anchor)
                .cloned()
                .unwrap_or_else(|| anchor.to_string()),
            anchor_rank_code: self.node(anchor).map(|n| n.rank_code.clone()),
            n_taxa,
            n_leaves,
            max_depth,
            mean_depth,
            topology_scale,
            max_children,
            branching_taxa_fraction,
            top_child_fraction,
        })
    }

    pub fn topology_stat(
        &self,
        taxon: TaxonId,
        anchor_rank: Option<&str>,
        stat: TopologyStat,
    ) -> Result<f64> {
        let p = self.topology(taxon, anchor_rank)?;
        Ok(match stat {
            TopologyStat::NTaxa => p.n_taxa as f64,
            TopologyStat::NLeaves => p.n_leaves as f64,
            TopologyStat::MaxDepth => p.max_depth as f64,
            TopologyStat::MeanDepth => p.mean_depth,
            TopologyStat::TopologyScale => p.topology_scale as f64,
            TopologyStat::MaxChildren => p.max_children as f64,
            TopologyStat::BranchingTaxaFraction => p.branching_taxa_fraction,
            TopologyStat::TopChildFraction => p.top_child_fraction,
        })
    }

    /// Batch form of [`Self::topology_stat`], evaluated in parallel in input order.
    pub fn topology_stats(
        &self,
        taxa: &[TaxonId],
        anchor_rank: Option<&str>,
        stat: TopologyStat,
    ) -> Result<Vec<f64>> {
        taxa.par_iter()
            .map(|taxon| self.topology_stat(*taxon, anchor_rank, stat))
            .collect()
    }
}

fn build_indexes(
    parent: &HashMap<TaxonId, Option<TaxonId>>,
    children: &HashMap<TaxonId, Vec<TaxonId>>,
) -> (HashMap<TaxonId, usize>, HashMap<TaxonId, (usize, usize)>) {
    let mut depth = HashMap::new();
    let mut intervals = HashMap::new();
    let mut visited = HashSet::new();
    let mut tick: usize = 0;
    let mut roots = parent
        .iter()
        .filter_map(|(node, p)| {
            if p.is_none()
                || *p == Some(*node)
                || p.is_some_and(|value| !parent.contains_key(&value))
            {
                Some(*node)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    roots.sort_unstable();
    roots.extend(parent.keys().copied());
    for root in roots {
        if visited.contains(&root) {
            continue;
        }
        let mut stack = vec![(root, 0, false)];
        while let Some((node, node_depth, exiting)) = stack.pop() {
            if exiting {
                if let Some((start, _)) = intervals.get(&node).copied() {
                    intervals.insert(node, (start, tick.saturating_sub(1)));
                }
                continue;
            }
            if !visited.insert(node) {
                continue;
            }
            intervals.insert(node, (tick, tick));
            depth.insert(node, node_depth);
            tick += 1;
            stack.push((node, node_depth, true));
            if let Some(values) = children.get(&node) {
                for child in values.iter().rev() {
                    stack.push((*child, node_depth + 1, false));
                }
            }
        }
    }
    (depth, intervals)
}

pub(crate) fn rank_to_code(rank: &str) -> Result<char> {
    let rank = rank.trim().to_ascii_uppercase();
    match rank.as_str() {
        "U" | "UNCLASSIFIED" => Ok('U'),
        "R" | "ROOT" => Ok('R'),
        "D" | "DOMAIN" | "SUPERKINGDOM" | "REALM" => Ok('D'),
        "K" | "KINGDOM" => Ok('K'),
        "P" | "PHYLUM" => Ok('P'),
        "C" | "CLASS" | "CLADE" => Ok('C'),
        "O" | "ORDER" => Ok('O'),
        "F" | "FAMILY" | "SUBFAMILY" => Ok('F'),
        "G" | "GENUS" => Ok('G'),
        "S" | "SPECIES" => Ok('S'),
        _ => bail!("rank must be one of: {}", RANK_CODES.join(", ")),
    }
}

pub(crate) fn rank_index(code: char) -> u8 {
    RANK_CODES
        .iter()
        .position(|candidate| candidate.starts_with(code))
        .unwrap_or(0) as u8
}

pub(crate) fn canonical_name(code: char) -> &'static str {
    match code {
        'U' => "unclassified",
        'R' => "root",
        'D' => "domain",
        'K' => "kingdom",
        'P' => "phylum",
        'C' => "class",
        'O' => "order",
        'F' => "family",
        'G' => "genus",
        'S' => "species",
        _ => "unclassified",
    }
}

pub(crate) fn major_rank_code(rank: &str) -> Option<char> {
    match rank {
        "root" | "acellular root" | "cellular root" => Some('R'),
        "domain" | "superkingdom" | "realm" => Some('D'),
        "kingdom" => Some('K'),
        "phylum" => Some('P'),
        "class" => Some('C'),
        "order" => Some('O'),
        "family" => Some('F'),
        "genus" => Some('G'),
        "species" => Some('S'),
        _ => None,
    }
}

pub(crate) fn assign_rank_codes(
    parent: &HashMap<TaxonId, Option<TaxonId>>,
    ranks: &HashMap<TaxonId, String>,
) -> HashMap<TaxonId, String> {
    fn one(
        taxon: TaxonId,
        parent: &HashMap<TaxonId, Option<TaxonId>>,
        ranks: &HashMap<TaxonId, String>,
        cache: &mut HashMap<TaxonId, String>,
        visiting: &mut HashSet<TaxonId>,
    ) -> String {
        if let Some(value) = cache.get(&taxon) {
            return value.clone();
        }
        if !visiting.insert(taxon) {
            cache.insert(taxon, "R".to_owned());
            return "R".to_owned();
        }
        let code = if taxon == 0 {
            "U".to_owned()
        } else if taxon == 1 {
            "R".to_owned()
        } else {
            let raw = ranks.get(&taxon).and_then(|rank| major_rank_code(rank));
            match parent.get(&taxon).copied().flatten() {
                None => raw.unwrap_or('R').to_string(),
                Some(p) if p == taxon || !parent.contains_key(&p) => raw.unwrap_or('R').to_string(),
                Some(p) => {
                    let parent_code = one(p, parent, ranks, cache, visiting);
                    let parent_base = parent_code.chars().next().unwrap_or('R');
                    if raw.is_some_and(|r| rank_index(r) > rank_index(parent_base)) {
                        raw.unwrap().to_string()
                    } else {
                        let depth = parent_code
                            .get(1..)
                            .and_then(|s| s.parse::<u32>().ok())
                            .unwrap_or(1);
                        format!("{}{}", parent_base, depth + 1)
                    }
                }
            }
        };
        visiting.remove(&taxon);
        cache.insert(taxon, code.clone());
        code
    }
    let mut cache = HashMap::new();
    for taxon in parent.keys() {
        one(*taxon, parent, ranks, &mut cache, &mut HashSet::new());
    }
    cache
}

pub(crate) fn taxonomic_order(
    present: HashSet<TaxonId>,
    parent: &HashMap<TaxonId, Option<TaxonId>>,
    ranks: &HashMap<TaxonId, String>,
    names: &HashMap<TaxonId, String>,
) -> Vec<TaxonId> {
    let mut ancestors = HashSet::new();
    let mut stack = present.iter().copied().collect::<Vec<_>>();
    while let Some(taxon) = stack.pop() {
        if let Some(p) = parent.get(&taxon).copied().flatten()
            && ancestors.insert(p)
        {
            stack.push(p);
        }
    }
    let nodes = present.union(&ancestors).copied().collect::<HashSet<_>>();
    let mut children: HashMap<TaxonId, Vec<TaxonId>> = HashMap::new();
    for taxon in &nodes {
        if let Some(p) = parent
            .get(taxon)
            .copied()
            .flatten()
            .filter(|p| nodes.contains(p))
        {
            children.entry(p).or_default().push(*taxon);
        }
    }
    let key = |taxon: &TaxonId| {
        (
            ranks.get(taxon).cloned().unwrap_or_default(),
            names.get(taxon).cloned().unwrap_or_default(),
            *taxon,
        )
    };
    for values in children.values_mut() {
        values.sort_by_key(&key);
    }
    let special = [0, 1, 9606, 2, 10239];
    let mut roots = nodes
        .iter()
        .filter(|taxon| {
            parent
                .get(taxon)
                .copied()
                .flatten()
                .is_none_or(|p| !nodes.contains(&p))
        })
        .copied()
        .collect::<Vec<_>>();
    roots.sort_by_key(|taxon| {
        (
            special
                .iter()
                .position(|v| v == taxon)
                .unwrap_or(usize::MAX),
            key(taxon),
        )
    });
    let mut order = Vec::new();
    let mut seen = HashSet::new();
    for root in roots {
        let mut dfs = vec![(root, false)];
        while let Some((node, expanded)) = dfs.pop() {
            if expanded {
                continue;
            }
            if !seen.insert(node) {
                continue;
            }
            if present.contains(&node) {
                order.push(node);
            }
            if let Some(values) = children.get(&node) {
                dfs.extend(values.iter().rev().map(|child| (*child, false)));
            }
        }
    }
    for taxon in BTreeSet::from_iter(present) {
        if !seen.contains(&taxon) {
            order.push(taxon);
        }
    }
    order
}

#[cfg(test)]
mod tests {
    use super::*;

    pub(crate) fn fixture() -> TaxonomicUtils {
        let raw = [
            (1, None, "root", "root"),
            (2, Some(1), "superkingdom", "Bacteria"),
            (10, Some(2), "family", "Family"),
            (11, Some(10), "genus", "Genus"),
            (12, Some(11), "species", "Species A"),
            (13, Some(11), "species", "Species B"),
        ];
        let parent = raw.iter().map(|(t, p, _, _)| (*t, *p)).collect();
        let ranks = raw
            .iter()
            .map(|(t, _, r, _)| (*t, (*r).to_owned()))
            .collect();
        let codes = assign_rank_codes(&parent, &ranks);
        let nodes = raw
            .iter()
            .map(|(taxon, parent, rank, _)| {
                let rank_code = codes[taxon].clone();
                let rank_base = rank_code.chars().next().unwrap();
                TaxonNode {
                    taxon: *taxon,
                    parent: *parent,
                    rank: (*rank).to_owned(),
                    rank_code,
                    rank_base,
                    rank_idx: rank_index(rank_base),
                    new_rank: canonical_name(rank_base).to_owned(),
                }
            })
            .collect();
        TaxonomicUtils::from_parts(
            raw.iter()
                .map(|(t, _, _, n)| (*t, (*n).to_owned()))
                .collect(),
            nodes,
            vec![],
            HashMap::new(),
            true,
            false,
            ".".into(),
        )
    }

    #[test]
    fn tree_queries_match_reference() {
        let tu = fixture();
        assert_eq!(tu.get_branch(12), vec![1, 2, 10, 11, 12]);
        assert_eq!(tu.get_subtree(11), vec![11, 12, 13]);
        assert!(tu.is_child(12, 11));
        assert!(tu.is_descendent(12, 10));
        assert!(!tu.is_descendent(10, 10));
        assert_eq!(tu.are_leaves(&[11, 12, 13]), vec![false, true, true]);
        assert_eq!(
            tu.are_children(&[12, 13], &[11, 11]).unwrap(),
            vec![true, true]
        );
        assert_eq!(
            tu.are_descendents(&[12, 10], &[10, 10]).unwrap(),
            vec![true, false]
        );
        assert!(tu.are_children(&[12], &[11, 10]).is_err());
        assert_eq!(tu.get_lca(12, 13), 11);
        assert_eq!(tu.get_distance(12, 13), 2);
        assert_eq!(tu.get_ancestor(12, "F").unwrap(), 10);
    }

    #[test]
    fn topology_matches_reference_formulae() {
        let p = fixture().topology(11, None).unwrap();
        assert_eq!(p.n_taxa, 3);
        assert_eq!(p.n_leaves, 2);
        assert_eq!(p.max_depth, 1);
        assert_eq!(p.mean_depth, 2.0 / 3.0);
        assert_eq!(p.top_child_fraction, 0.5);
    }
}
