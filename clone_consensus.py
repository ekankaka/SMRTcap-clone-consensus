#!/usr/bin/env python3

import argparse
import csv
import json
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

VERSION = "1.0.0"
REQUIRED = ["participant_id", "sample_id", "read", "clone_id", "STRAND", "HIV_SEQ"]
VALID_DNA = set("ACGTNRYSWKMBDHV")
BAD_ID = re.compile(r"[\s,;|]")


def clean(value):
    return str(value).strip()


def normalize(seq, strand):
    seq = "".join(str(seq).split()).upper()
    if not seq:
        raise ValueError("HIV_SEQ cannot be blank")
    bad = set(seq) - VALID_DNA
    if bad:
        raise ValueError("HIV_SEQ contains unsupported character(s): " + "".join(sorted(bad)))
    if strand == "plus":
        return seq
    if strand == "minus":
        return str(Seq(seq).reverse_complement())
    raise ValueError(f"STRAND must be plus or minus, not: {strand}")


def one_reference(path):
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        raise ValueError("Reference FASTA must contain exactly one sequence")
    seq = str(records[0].seq).upper().replace("-", "")
    if not seq:
        raise ValueError("Reference sequence is empty")
    bad = set(seq) - VALID_DNA
    if bad:
        raise ValueError("Reference contains unsupported character(s): " + "".join(sorted(bad)))
    return seq


def prepare(args):
    work = Path(args.work)
    clones_dir = work / "clones"
    work.mkdir(parents=True, exist_ok=True)
    clones_dir.mkdir(parents=True, exist_ok=True)

    ref = one_reference(args.reference)
    SeqIO.write([SeqRecord(Seq(ref), id="REFERENCE", description="")], work / "reference.fasta", "fasta")

    with open(args.input, newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError("Input CSV has no header")
        missing = [x for x in REQUIRED if x not in reader.fieldnames]
        if missing:
            raise ValueError("Missing required column(s): " + ", ".join(missing))
        rows = list(reader)
    if not rows:
        raise ValueError("Input CSV contains no data rows")

    normalized_rows = []
    clones = defaultdict(list)
    seen = set()

    for line_no, row in enumerate(rows, start=2):
        pid = clean(row["participant_id"])
        sid = clean(row["sample_id"])
        read = clean(row["read"])
        cid = clean(row["clone_id"])
        strand = clean(row["STRAND"]).lower()

        for name, value in [("participant_id", pid), ("sample_id", sid), ("read", read), ("clone_id", cid)]:
            if not value:
                raise ValueError(f"Blank {name} at CSV line {line_no}")
            if BAD_ID.search(value):
                raise ValueError(f"{name} cannot contain whitespace, comma, semicolon, or '|': {value}")

        sequence_id = f"{pid}|{sid}|{read}"
        if sequence_id in seen:
            raise ValueError(f"Duplicate participant_id + sample_id + read: {sequence_id}")
        seen.add(sequence_id)

        raw = "".join(str(row["HIV_SEQ"]).split()).upper()
        forward = normalize(raw, strand)
        out = {
            "participant_id": pid, "sample_id": sid, "read": read, "clone_id": cid,
            "STRAND": strand, "HIV_SEQ": raw, "HIV_SEQ_FORWARD": forward,
            "sequence_id": sequence_id,
        }
        normalized_rows.append(out)
        clones[(pid, cid)].append(out)

    fields = REQUIRED + ["HIV_SEQ_FORWARD", "sequence_id"]
    with open(work / "normalized_sequences.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(normalized_rows)

    manifest = []
    for i, ((pid, cid), items) in enumerate(sorted(clones.items()), start=1):
        key = f"clone_{i:06d}"
        folder = clones_dir / key
        folder.mkdir(parents=True, exist_ok=True)
        SeqIO.write([
            SeqRecord(Seq(x["HIV_SEQ_FORWARD"]), id=x["sequence_id"], description="") for x in items
        ], folder / "fragments.fasta", "fasta")

        sample_ids = sorted({x["sample_id"] for x in items})
        lengths = [len(x["HIV_SEQ_FORWARD"]) for x in items]
        info = {
            "clone_key": key, "participant_id": pid, "clone_id": cid,
            "n_fragments": len(items), "n_samples": len(sample_ids), "sample_ids": sample_ids,
            "shortest_fragment_nt": min(lengths), "longest_fragment_nt": max(lengths),
        }
        with open(folder / "clone_info.json", "w") as handle:
            json.dump(info, handle, indent=2)
        manifest.append({**info, "sample_ids": ";".join(sample_ids)})

    manifest_fields = [
        "clone_key", "participant_id", "clone_id", "n_fragments", "n_samples", "sample_ids",
        "shortest_fragment_nt", "longest_fragment_nt",
    ]
    with open(work / "clone_manifest.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields)
        writer.writeheader()
        writer.writerows(manifest)

    print(f"Prepared {len(normalized_rows)} fragments in {len(manifest)} clones")


def interval(seq):
    positions = [i for i, base in enumerate(seq) if base != "-"]
    return (positions[0], positions[-1]) if positions else None


def consensus(args):
    alignment = Path(args.alignment)
    folder = alignment.parent
    records = list(SeqIO.parse(alignment, "fasta"))
    refs = [r for r in records if r.id == "REFERENCE"]
    fragments = [r for r in records if r.id != "REFERENCE"]
    if len(refs) != 1 or not fragments:
        raise ValueError("Alignment must contain one REFERENCE and at least one fragment")

    length = len(refs[0].seq)
    if any(len(r.seq) != length for r in records):
        raise ValueError("Alignment sequences have different lengths")

    seqs = [str(r.seq).upper() for r in fragments]
    intervals = [interval(s) for s in seqs]
    calls, depths = [], []
    disagreements = 0

    for col in range(length):
        votes = []
        for seq, span in zip(seqs, intervals):
            if span and span[0] <= col <= span[1] and seq[col] in "ACGT-":
                votes.append(seq[col])
        counts = Counter(votes)
        depth = len(votes)
        depths.append(depth)
        if len(counts) > 1:
            disagreements += 1
        if depth == 0:
            calls.append("N")
        else:
            base, count = counts.most_common(1)[0]
            calls.append(base if count > depth / 2 else "N")

    covered = [i for i, d in enumerate(depths) if d > 0]
    if not covered:
        raise ValueError("No informative fragment coverage in alignment")
    first, last = covered[0], covered[-1]
    sequence = "".join(x for x in calls[first:last + 1] if x != "-")
    if not sequence:
        raise ValueError("Consensus sequence is empty")

    with open(folder / "clone_info.json") as handle:
        info = json.load(handle)
    consensus_id = f"{info['participant_id']}|{info['clone_id']}|n={info['n_fragments']}"
    SeqIO.write([SeqRecord(Seq(sequence), id=consensus_id, description="")], folder / "consensus.fasta", "fasta")

    positive_depth = [d for d in depths[first:last + 1] if d > 0]
    stats = {
        **info,
        "consensus_id": consensus_id,
        "consensus_type": "single_fragment" if info["n_fragments"] == 1 else "multi_fragment",
        "consensus_length_nt": len(sequence),
        "consensus_non_N_nt": sum(b in "ACGT" for b in sequence),
        "percent_consensus_N": round(100 * sequence.count("N") / len(sequence), 2),
        "mean_depth": round(statistics.mean(positive_depth), 2),
        "median_depth": round(statistics.median(positive_depth), 2),
        "min_depth": min(positive_depth),
        "max_depth": max(positive_depth),
        "disagreement_positions": disagreements,
    }
    with open(folder / "consensus_stats.json", "w") as handle:
        json.dump(stats, handle, indent=2)

    reference = str(refs[0].seq).upper()
    ref_pos = sum(1 for b in reference[:first] if b != "-")
    cons_pos = 0
    depth_fields = ["alignment_column", "reference_position", "consensus_position", "depth", "consensus_call"]
    with open(folder / "depth.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=depth_fields)
        writer.writeheader()
        for col in range(first, last + 1):
            if reference[col] != "-":
                ref_pos += 1
            call = calls[col]
            if call != "-":
                cons_pos += 1
            writer.writerow({
                "alignment_column": col + 1,
                "reference_position": ref_pos if reference[col] != "-" else "",
                "consensus_position": cons_pos if call != "-" else "",
                "depth": depths[col],
                "consensus_call": call,
            })


def best_paf(path):
    best = {}
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return best
    with open(path) as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            if len(f) < 12:
                continue
            qlen, qstart, qend = int(f[1]), int(f[2]), int(f[3])
            nmatch, aln_len, mapq = int(f[9]), int(f[10]), int(f[11])
            hit = {
                "strand": f[4],
                "alignment_bp": aln_len,
                "identity": nmatch / aln_len if aln_len else 0,
                "query_coverage": (qend - qstart) / qlen if qlen else 0,
                "mapq": mapq,
            }
            score = (hit["query_coverage"], hit["identity"], aln_len, mapq)
            if f[0] not in best or score > best[f[0]][0]:
                best[f[0]] = (score, hit)
    return {q: x[1] for q, x in best.items()}


def finalize(args):
    work, results = Path(args.work), Path(args.results)
    results.mkdir(parents=True, exist_ok=True)
    with open(work / "normalized_sequences.csv", newline="") as handle:
        normalized = list(csv.DictReader(handle))

    clone_dirs = sorted((work / "clones").glob("clone_*"))
    if not clone_dirs:
        raise ValueError("No clone work folders found")

    consensus_records, summaries, hits = [], [], {}
    for folder in clone_dirs:
        consensus_records.extend(SeqIO.parse(folder / "consensus.fasta", "fasta"))
        with open(folder / "consensus_stats.json") as handle:
            summaries.append(json.load(handle))
        hits.update(best_paf(folder / "fragments_vs_consensus.paf"))
    SeqIO.write(consensus_records, results / "clone_consensus.fasta", "fasta")

    fragment_rows = []
    clone_hits = defaultdict(list)
    for row in normalized:
        hit = hits.get(row["sequence_id"])
        fragment_rows.append({
            "participant_id": row["participant_id"],
            "sample_id": row["sample_id"],
            "read": row["read"],
            "clone_id": row["clone_id"],
            "sequence_id": row["sequence_id"],
            "fragment_length_nt": len(row["HIV_SEQ_FORWARD"]),
            "mapped": "yes" if hit else "no",
            "strand_to_consensus": hit["strand"] if hit else "",
            "alignment_bp": hit["alignment_bp"] if hit else "",
            "percent_identity": round(100 * hit["identity"], 2) if hit else "",
            "percent_fragment_covered": round(100 * hit["query_coverage"], 2) if hit else "",
            "mapq": hit["mapq"] if hit else "",
        })
        clone_hits[(row["participant_id"], row["clone_id"])].append(hit)

    qc_fields = [
        "participant_id", "sample_id", "read", "clone_id", "sequence_id", "fragment_length_nt",
        "mapped", "strand_to_consensus", "alignment_bp", "percent_identity", "percent_fragment_covered", "mapq",
    ]
    with open(results / "fragment_qc.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=qc_fields)
        writer.writeheader()
        writer.writerows(fragment_rows)

    for stats in summaries:
        valid_hits = [h for h in clone_hits[(stats["participant_id"], stats["clone_id"])] if h]
        stats["n_fragments_mapped_to_consensus"] = len(valid_hits)
        if valid_hits:
            identities = [100 * h["identity"] for h in valid_hits]
            coverages = [100 * h["query_coverage"] for h in valid_hits]
            stats["mean_remap_identity"] = round(statistics.mean(identities), 2)
            stats["min_remap_identity"] = round(min(identities), 2)
            stats["mean_fragment_coverage"] = round(statistics.mean(coverages), 2)
            stats["min_fragment_coverage"] = round(min(coverages), 2)
        else:
            for name in [
                "mean_remap_identity", "min_remap_identity", "mean_fragment_coverage", "min_fragment_coverage"
            ]:
                stats[name] = ""
        stats["sample_ids"] = ";".join(stats["sample_ids"])

    summary_fields = [
        "participant_id", "clone_id", "consensus_id", "consensus_type", "n_fragments", "n_samples",
        "sample_ids", "shortest_fragment_nt", "longest_fragment_nt", "consensus_length_nt",
        "consensus_non_N_nt", "percent_consensus_N", "mean_depth", "median_depth", "min_depth", "max_depth",
        "disagreement_positions", "n_fragments_mapped_to_consensus", "mean_remap_identity", "min_remap_identity",
        "mean_fragment_coverage", "min_fragment_coverage",
    ]
    with open(results / "clone_consensus_summary.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summaries)

    print(f"Consensus sequences: {results / 'clone_consensus.fasta'}")
    print(f"Clone summary:       {results / 'clone_consensus_summary.csv'}")
    print(f"Fragment QC:         {results / 'fragment_qc.csv'}")


def parser():
    p = argparse.ArgumentParser(description="SMRTcap clone consensus helper")
    sub = p.add_subparsers(dest="command", required=True)

    a = sub.add_parser("prepare", help="validate input, normalize strands, and create clone FASTAs")
    a.add_argument("--input", required=True)
    a.add_argument("--reference", required=True)
    a.add_argument("--work", required=True)
    a.set_defaults(func=prepare)

    a = sub.add_parser("consensus", help="call a coverage-aware consensus from one clone alignment")
    a.add_argument("--alignment", required=True)
    a.set_defaults(func=consensus)

    a = sub.add_parser("finalize", help="combine clone consensuses and summarize remapping QC")
    a.add_argument("--work", required=True)
    a.add_argument("--results", required=True)
    a.set_defaults(func=finalize)

    a = sub.add_parser("version", help="print pipeline version")
    a.set_defaults(func=lambda args: print(VERSION))
    return p


if __name__ == "__main__":
    args = parser().parse_args()
    args.func(args)
