from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


TM_SCORES = Path(
    "/Users/dalsasso/Desktop/Posdoc/CAU/network_analysis/Zpa796/tm-score_network/"
    "pairwise_TM0.5/Zpa796_tm-scores_all-vs-all_summary.txt"
)
ANNOTATION = Path(
    "/Users/dalsasso/Desktop/Posdoc/CAU/annotations/table-annot_backup/"
    "Zpa796_secretome-annotation_12-08-24.tsv"
)
ORIGINAL_GRAPH = Path(
    "/Users/dalsasso/Desktop/Posdoc/CAU/network_analysis/Zpa796/tm-score_network/"
    "pairwise_TM0.5/Zpa796_tm-scores_network.graphml"
)
OUTDIR = Path(
    "/Users/dalsasso/Desktop/Posdoc/CAU/network_analysis/Zpa796/tm-score_network/"
    "pairwise_TM0.5/pLDDT_threshold_sensitivity_results"
)


TM_THRESHOLD = 0.5
THRESHOLDS = list(range(50, 100, 5))
FIGURE_SIZE = (8, 6)
LINE_AXIS_LABEL_SIZE = 18
LINE_TICK_SIZE = 16
LINE_LEGEND_SIZE = 16


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def map_effector_status(value):
    if pd.isna(value):
        return np.nan
    if "Non-effector" in str(value):
        return "Non-effector"
    return "Effector"


def read_qualifying_edges(path):
    pair_scores = {}
    with path.open() as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            protein1, protein2, score_text = parts
            try:
                score = float(score_text)
            except ValueError:
                continue
            pair = tuple(sorted((protein1, protein2)))
            pair_scores.setdefault(pair, []).append(score)

    return {
        pair: max(scores)
        for pair, scores in pair_scores.items()
        if len(scores) == 2 and all(score >= TM_THRESHOLD for score in scores)
    }


def build_network(edges, allowed_proteins):
    graph = nx.Graph()
    for (protein1, protein2), score in edges.items():
        if protein1 in allowed_proteins and protein2 in allowed_proteins:
            graph.add_edge(protein1, protein2, weight=score)
    return graph


def sorted_components(graph):
    return sorted(
        nx.connected_components(graph),
        key=lambda component: (-len(component), sorted(component)[0]),
    )


def component_assignments(graph, annotation, rebuilt=True):
    rows = []
    if rebuilt:
        for index, component in enumerate(sorted_components(graph), start=1):
            rebuilt_label = f"G.{index:02d}"
            parents = sorted(
                annotation.loc[list(component), "Structural subgraph"].dropna().unique()
            )
            parent_label = ",".join(parents)
            for protein in sorted(component):
                rows.append(
                    {
                        "Protein ID": protein,
                        "Rebuilt structural subgraph": rebuilt_label,
                        "Parent structural subgraph": parent_label,
                        "Rebuilt cluster size": len(component),
                    }
                )
    else:
        for protein in graph.nodes():
            label = annotation.at[protein, "Structural subgraph"]
            if pd.notna(label):
                rows.append(
                    {
                        "Protein ID": protein,
                        "Rebuilt structural subgraph": label,
                        "Parent structural subgraph": label,
                    }
                )
    return pd.DataFrame(rows)


def enrichment_input(annotation, assignments):
    return annotation.reset_index(drop=True).merge(assignments, on="Protein ID", how="inner")


def run_enrichment(data, group_column, status_column, positive_label):
    filtered = data[[group_column, status_column]].dropna()
    contingency = pd.crosstab(filtered[group_column], filtered[status_column])
    negative_labels = [column for column in contingency.columns if column != positive_label]
    if positive_label not in contingency.columns or len(negative_labels) != 1:
        raise ValueError(
            f"Unexpected categories for {status_column}: {contingency.columns.tolist()}"
        )
    negative_label = negative_labels[0]
    contingency = contingency[[positive_label, negative_label]]

    total_positive = int(contingency[positive_label].sum())
    total_negative = int(contingency[negative_label].sum())
    rows = []
    for cluster, row in contingency.iterrows():
        positive = int(row[positive_label])
        negative = int(row[negative_label])
        table = [
            [positive, negative],
            [total_positive - positive, total_negative - negative],
        ]
        odds_ratio, p_value = fisher_exact(table, alternative="greater")
        rows.append(
            {
                "Rebuilt structural subgraph": cluster,
                "Odds Ratio": odds_ratio,
                "P-value": p_value,
                "Cluster Size": positive + negative,
                "Positive in cluster": positive,
                "Negative in cluster": negative,
                "Total positives": total_positive,
                "Total negatives": total_negative,
            }
        )

    results = pd.DataFrame(rows)
    results["Adjusted P-value"] = multipletests(results["P-value"], method="fdr_bh")[1]
    results["Significant"] = results["Adjusted P-value"] < 0.05
    return results.sort_values(["Adjusted P-value", "P-value"])


def run_enrichments(annotation, assignments):
    data = enrichment_input(annotation, assignments)
    analyses = {
        "AMP": ("AM prediction (AF2/AMAPEC)", "Antimicrobial"),
        "Effector": ("Effector status", "Effector"),
    }
    results = {}
    parent_map = assignments.drop_duplicates("Rebuilt structural subgraph").set_index(
        "Rebuilt structural subgraph"
    )["Parent structural subgraph"]
    for analysis, (column, positive_label) in analyses.items():
        result = run_enrichment(
            data, "Rebuilt structural subgraph", column, positive_label
        )
        result["Parent structural subgraph"] = result[
            "Rebuilt structural subgraph"
        ].map(parent_map)
        results[analysis] = result
    return results


def original_effector_enriched_groups(original_graph, annotation):
    assignments = component_assignments(original_graph, annotation, rebuilt=False)
    results = run_enrichments(annotation, assignments)["Effector"]
    results.to_csv(OUTDIR / "original_effector_enrichment.tsv", sep="\t", index=False)
    return sorted(
        results.loc[results["Significant"], "Parent structural subgraph"].unique()
    )


def original_members_by_group(original_groups, original_graph, annotation):
    original_nodes = set(original_graph.nodes())
    members = {}
    for group in original_groups:
        group_members = set(annotation.index[annotation["Structural subgraph"] == group])
        members[group] = group_members & original_nodes
    return members


def cluster_jaccard_rows(original_members, graph, threshold):
    rebuilt_components = list(nx.connected_components(graph))
    graph_nodes = set(graph.nodes())
    rows = []
    for group, members in original_members.items():
        best_jaccard = 0
        best_size = 0
        retained = members & graph_nodes
        for component in rebuilt_components:
            component = set(component)
            intersection = members & component
            if intersection:
                jaccard = len(intersection) / len(members | component)
                if jaccard > best_jaccard:
                    best_jaccard = jaccard
                    best_size = len(intersection)
        rows.append(
            {
                "Threshold": threshold,
                "Original effector-enriched group": group,
                "Original cluster size": len(members),
                "Retained proteins": len(retained),
                "Member retention proportion": len(retained) / len(members),
                "Best-component overlap": best_size,
                "Jaccard similarity": best_jaccard,
            }
        )
    return rows


def plot_normalized_network_retention(summary, original_graph, output_stem):
    fig, axis = plt.subplots(figsize=FIGURE_SIZE)
    normalized = {
        "Nodes": summary["Nodes"] / original_graph.number_of_nodes(),
        "Edges": summary["Edges"] / original_graph.number_of_edges(),
        "Structural clusters": summary["Connected components"]
        / nx.number_connected_components(original_graph),
    }
    markers = ["o", "s", "^"]
    for (label, values), color, marker in zip(
        normalized.items(), ["#08306B", "#2171B5", "#6BAED6"], markers
    ):
        axis.plot(
            summary["Threshold"],
            values,
            marker=marker,
            markersize=7,
            color=color,
            label=label,
        )
    axis.set_xlabel("Minimum mean pLDDT (AF2)", fontsize=LINE_AXIS_LABEL_SIZE)
    axis.set_ylabel("Proportion retained", fontsize=LINE_AXIS_LABEL_SIZE)
    axis.set_ylim(0, 1.05)
    axis.set_xticks(THRESHOLDS)
    axis.tick_params(axis="both", labelsize=LINE_TICK_SIZE)
    axis.legend(frameon=False, fontsize=LINE_LEGEND_SIZE, loc="lower left")
    axis.grid(axis="y", alpha=0.25)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_stem.with_suffix(".pdf"))
    plt.close(fig)


def plot_enriched_clusters(summary, output_stem):
    fig, axis = plt.subplots(figsize=FIGURE_SIZE)
    axis.plot(
        summary["Threshold"],
        summary["Original effector groups remaining significant"],
        marker="o",
        markersize=7,
        color="#C55A11",
        label="Effector-enriched clusters",
    )
    axis.set_xlabel("Minimum mean pLDDT (AF2)", fontsize=LINE_AXIS_LABEL_SIZE)
    axis.set_ylabel("Clusters retained", fontsize=LINE_AXIS_LABEL_SIZE)
    axis.set_xticks(THRESHOLDS)
    axis.set_ylim(bottom=0)
    axis.tick_params(axis="both", labelsize=LINE_TICK_SIZE)
    axis.legend(frameon=False, fontsize=LINE_LEGEND_SIZE, loc="lower left")
    axis.grid(axis="y", alpha=0.25)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_stem.with_suffix(".pdf"))
    plt.close(fig)


def plot_member_retention_heatmap(stability, significant, output_stem):
    matrix = stability.pivot(
        index="Original effector-enriched group",
        columns="Threshold",
        values="Member retention proportion",
    ).sort_index()
    fig, axis = plt.subplots(figsize=FIGURE_SIZE)
    axis = sns.heatmap(
        matrix,
        annot=False,
        cmap="Blues",
        cbar=False,
        vmin=0,
        vmax=1,
        square=True,
        ax=axis,
    )
    significant_effector = significant.loc[
        significant["Analysis"].eq("Effector")
        & significant["Parent structural subgraph"].isin(matrix.index),
        ["Threshold", "Parent structural subgraph"],
    ]
    significant_cells = set(
        zip(
            significant_effector["Parent structural subgraph"],
            significant_effector["Threshold"],
        )
    )
    for row_index, group in enumerate(matrix.index):
        for column_index, threshold in enumerate(matrix.columns):
            if (group, threshold) in significant_cells:
                axis.scatter(
                    column_index + 0.5,
                    row_index + 0.5,
                    marker="D",
                    s=55,
                    facecolors="white",
                    edgecolors="black",
                    linewidths=0.5,
                    zorder=5,
                )
    axis.set_title("")
    axis.set_xlabel("Minimum mean pLDDT (AF2)", fontsize=22)
    axis.set_ylabel("")
    axis.set_xticklabels(axis.get_xticklabels(), rotation=0, ha="center", fontsize=18)
    axis.set_yticklabels(axis.get_yticklabels(), rotation=0, fontsize=18)
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)
        spine.set_edgecolor("black")
    fig.subplots_adjust(left=0.10, right=0.82, bottom=0.18, top=0.95)
    fig.canvas.draw()
    position = axis.get_position()
    colorbar_axis = fig.add_axes(
        [position.x1 + 0.025, position.y0, 0.025, position.height]
    )
    colorbar = fig.colorbar(axis.collections[0], cax=colorbar_axis)
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(labelsize=18)
    fig.savefig(output_stem.with_suffix(".pdf"))
    plt.close(fig)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    annotation = pd.read_csv(ANNOTATION, sep="\t").set_index("Protein ID", drop=False)
    annotation["Effector status"] = annotation["EffectorP3"].apply(map_effector_status)
    original_graph = nx.read_graphml(ORIGINAL_GRAPH)
    qualifying_edges = read_qualifying_edges(TM_SCORES)
    enriched_groups = original_effector_enriched_groups(original_graph, annotation)
    original_members = original_members_by_group(enriched_groups, original_graph, annotation)

    summary_rows = []
    significant_rows = []
    stability_rows = []

    for threshold in THRESHOLDS:
        allowed = set(annotation.index[annotation["pLDDT (AF2)"] >= threshold])
        graph = build_network(qualifying_edges, allowed)
        assignments = component_assignments(graph, annotation)

        enrichments = run_enrichments(annotation, assignments)
        counts = {}
        original_effector_groups_remaining = 0
        for analysis, results in enrichments.items():
            counts[analysis] = int(results["Significant"].sum())
            if analysis == "Effector":
                original_effector_groups_remaining = results.loc[
                    results["Significant"]
                    & results["Parent structural subgraph"].isin(enriched_groups),
                    "Parent structural subgraph",
                ].nunique()
            significant = results[results["Significant"]].copy()
            significant.insert(0, "Analysis", analysis)
            significant.insert(0, "Threshold", threshold)
            significant_rows.extend(significant.to_dict("records"))

        summary_rows.append(
            {
                "Threshold": threshold,
                "Input proteins": len(allowed),
                "Nodes": graph.number_of_nodes(),
                "Edges": graph.number_of_edges(),
                "Connected components": nx.number_connected_components(graph),
                "Effector enriched clusters": counts["Effector"],
                "Original effector groups remaining significant": (
                    original_effector_groups_remaining
                ),
                "AMP enriched clusters": counts["AMP"],
            }
        )
        stability_rows.extend(cluster_jaccard_rows(original_members, graph, threshold))

    summary = pd.DataFrame(summary_rows)
    significant = pd.DataFrame(significant_rows)
    stability = pd.DataFrame(stability_rows)
    summary.to_csv(OUTDIR / "threshold_summary.tsv", sep="\t", index=False)
    significant.to_csv(OUTDIR / "significant_enrichments.tsv", sep="\t", index=False)
    stability.to_csv(
        OUTDIR / "original_effector_cluster_member_retention.tsv",
        sep="\t",
        index=False,
    )

    plot_normalized_network_retention(
        summary, original_graph, OUTDIR / "normalized_network_retention"
    )
    plot_enriched_clusters(summary, OUTDIR / "enriched_clusters_by_threshold")
    plot_member_retention_heatmap(
        stability,
        significant,
        OUTDIR / "original_effector_cluster_member_retention_heatmap",
    )

    print("Original effector-enriched clusters:", ", ".join(enriched_groups))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
