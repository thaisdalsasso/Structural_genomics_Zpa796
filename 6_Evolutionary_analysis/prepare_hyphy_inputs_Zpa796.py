#!/usr/bin/env python3

'''

Prepares codon-aware alignments from protein MSAs, retrieves the matching CDS records,
threads codons through the protein alignments, and writes trimmed codon-aware alignments ready for HyPhy.

'''

from __future__ import annotations

import csv
import re
import textwrap
from pathlib import Path


ROOT = Path("/path2/phylogeny_homologs")
SOURCE_ROOT = Path("/path2/structural_clusters_MSAs")
MATURE_CDS_FASTA = Path("/data/References/cds/Zpa796_mature-CDS.secretome.fasta")
OUTDIR = ROOT / "hyphy_analysis_Zpa796"

CLUSTERS = {
    "G.04": {
        "mature_protein": SOURCE_ROOT / "G.04/mafft-linsi/G.04_mature.fasta",
        "protein_alignment": SOURCE_ROOT / "G.04/mafft-linsi/G.04_mature_mafft-linsi.fasta",
        "trimmed_protein_alignment": SOURCE_ROOT / "G.04/mafft-linsi/G.04_mature_mafft-linsi_trimmed.fasta",
    },
    "G.12": {
        "mature_protein": SOURCE_ROOT / "G.12/mafft-linsi/G.12_mature.fasta",
        "protein_alignment": SOURCE_ROOT / "G.12/mafft-linsi/G.12_mature_mafft-linsi.fasta",
        "trimmed_protein_alignment": SOURCE_ROOT / "G.12/mafft-linsi/G.12_mature_mafft-linsi_trimmed.fasta",
    },
}

GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "G",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Correct the single GCG entry after keeping the table compact above.
GENETIC_CODE["GCG"] = "A"


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header = None
    chunks: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records[header] = "".join(chunks)
    return records


def write_fasta(path: Path, records: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for seq_id, sequence in records.items():
            handle.write(f">{seq_id}\n")
            handle.write("\n".join(textwrap.wrap(sequence, width=80)))
            handle.write("\n")


def clean_cds(sequence: str) -> str:
    cds = re.sub(r"[^A-Za-z]", "", sequence).upper().replace("U", "T")
    if len(cds) >= 3 and GENETIC_CODE.get(cds[-3:], "X") == "*":
        cds = cds[:-3]
    return cds


def translate(cds: str) -> str:
    return "".join(GENETIC_CODE.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))


def thread_cds(protein_alignment: dict[str, str], cds_records: dict[str, str]) -> dict[str, str]:
    codon_alignment: dict[str, str] = {}
    for seq_id, aa_alignment in protein_alignment.items():
        cds = clean_cds(cds_records[seq_id])
        if len(cds) % 3 != 0:
            raise ValueError(f"{seq_id}: mature CDS length is not divisible by 3")
        codon_index = 0
        codons: list[str] = []
        for aa in aa_alignment:
            if aa == "-":
                codons.append("---")
            else:
                codon = cds[codon_index:codon_index + 3]
                if len(codon) != 3:
                    raise ValueError(f"{seq_id}: protein alignment consumes more codons than available")
                codons.append(codon)
                codon_index += 3
        if codon_index != len(cds):
            raise ValueError(f"{seq_id}: mature CDS has unused bases after threading")
        codon_alignment[seq_id] = "".join(codons)
    return codon_alignment


def retained_columns(untrimmed: dict[str, str], trimmed: dict[str, str]) -> list[int]:
    if set(untrimmed) != set(trimmed):
        raise ValueError("Untrimmed and trimmed protein alignments have different sequence IDs")
    seq_ids = sorted(untrimmed)
    untrimmed_len = len(next(iter(untrimmed.values())))
    trimmed_len = len(next(iter(trimmed.values())))
    kept: list[int] = []
    cursor = 0
    for trimmed_index in range(trimmed_len):
        target = tuple(trimmed[seq_id][trimmed_index] for seq_id in seq_ids)
        while cursor < untrimmed_len:
            current = tuple(untrimmed[seq_id][cursor] for seq_id in seq_ids)
            if current == target:
                kept.append(cursor)
                cursor += 1
                break
            cursor += 1
        else:
            raise ValueError(f"Could not map trimmed column {trimmed_index + 1}")
    return kept


def trim_codon_alignment(codon_alignment: dict[str, str], kept_columns: list[int]) -> dict[str, str]:
    trimmed: dict[str, str] = {}
    for seq_id, sequence in codon_alignment.items():
        trimmed[seq_id] = "".join(sequence[column * 3:column * 3 + 3] for column in kept_columns)
    return trimmed


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    all_mature_cds = {seq_id: clean_cds(seq) for seq_id, seq in read_fasta(MATURE_CDS_FASTA).items()}
    summary_rows: list[dict[str, object]] = []

    for cluster, paths in CLUSTERS.items():
        cluster_dir = OUTDIR / cluster
        cluster_dir.mkdir(parents=True, exist_ok=True)

        mature_proteins = read_fasta(paths["mature_protein"])
        protein_alignment = read_fasta(paths["protein_alignment"])
        trimmed_protein_alignment = read_fasta(paths["trimmed_protein_alignment"])
        missing = sorted(set(mature_proteins) - set(all_mature_cds))
        if missing:
            raise ValueError(f"{cluster}: missing mature CDS for {', '.join(missing)}")

        mature_cds = {seq_id: all_mature_cds[seq_id] for seq_id in mature_proteins}
        qa_rows = []
        for seq_id, protein in mature_proteins.items():
            translated = translate(mature_cds[seq_id])
            qa_rows.append({
                "Protein ID": seq_id,
                "Mature protein length": len(protein),
                "Mature CDS length": len(mature_cds[seq_id]),
                "Mature CDS codons": len(mature_cds[seq_id]) // 3,
                "Translation matches mature protein": translated == protein.upper().replace("*", ""),
            })

        codon_alignment = thread_cds(protein_alignment, mature_cds)
        kept = retained_columns(protein_alignment, trimmed_protein_alignment)
        trimmed_codon_alignment = trim_codon_alignment(codon_alignment, kept)

        write_fasta(cluster_dir / f"{cluster}_Zpa796_mature_proteins.fasta", mature_proteins)
        write_fasta(cluster_dir / f"{cluster}_Zpa796_mature_cds.fasta", mature_cds)
        write_fasta(cluster_dir / f"{cluster}_Zpa796_mature_protein_mafft_linsi.fasta", protein_alignment)
        write_fasta(cluster_dir / f"{cluster}_Zpa796_mature_protein_trimmed.fasta", trimmed_protein_alignment)
        write_fasta(cluster_dir / f"{cluster}_Zpa796_mature_codon_alignment.fasta", codon_alignment)
        write_fasta(cluster_dir / f"{cluster}_Zpa796_mature_codon_alignment_trimmed.fasta", trimmed_codon_alignment)

        with (cluster_dir / f"{cluster}_Zpa796_mature_sequence_QA.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(qa_rows[0].keys()))
            writer.writeheader()
            writer.writerows(qa_rows)

        with (cluster_dir / f"{cluster}_Zpa796_retained_protein_columns.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Retained column in trimmed alignment", "Original protein alignment column"])
            for trimmed_index, original_index in enumerate(kept, start=1):
                writer.writerow([trimmed_index, original_index + 1])

        summary_rows.append({
            "Cluster": cluster,
            "Sequences": len(mature_proteins),
            "Untrimmed protein columns": len(next(iter(protein_alignment.values()))),
            "Trimmed protein columns": len(next(iter(trimmed_protein_alignment.values()))),
            "Trimmed codon length": len(next(iter(trimmed_codon_alignment.values()))),
            "All translations match": all(row["Translation matches mature protein"] for row in qa_rows),
        })
        print(f"{cluster}: prepared {len(mature_proteins)} CDS/protein sequences")

    with (OUTDIR / "Zpa796_hyphy_input_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)



if __name__ == "__main__":
    main()
