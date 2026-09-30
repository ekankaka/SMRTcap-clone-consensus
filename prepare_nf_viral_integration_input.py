#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path

REQUIRED_OUTPUT = ["participant_id", "sample_id", "read", "clone_id", "STRAND", "HIV_SEQ"]
BAD_ID = re.compile(r"[\s,;|]")


def clean(value):
    return str(value).strip()


def load_mapping(path, sample_ids):
    if path is None:
        return {sample_id: sample_id for sample_id in sample_ids}

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError("Sample mapping CSV has no header")
        missing = [x for x in ["sample_id", "participant_id"] if x not in reader.fieldnames]
        if missing:
            raise ValueError("Missing mapping column(s): " + ", ".join(missing))
        rows = list(reader)

    mapping = {}
    for line_no, row in enumerate(rows, start=2):
        sample_id = clean(row["sample_id"])
        participant_id = clean(row["participant_id"])
        if not sample_id or not participant_id:
            raise ValueError(f"Blank sample_id or participant_id at mapping line {line_no}")
        if BAD_ID.search(sample_id) or BAD_ID.search(participant_id):
            raise ValueError("sample_id and participant_id cannot contain whitespace, comma, semicolon, or '|'")
        if sample_id in mapping:
            raise ValueError(f"Duplicate sample_id in mapping CSV: {sample_id}")
        mapping[sample_id] = participant_id

    missing_samples = sorted(set(sample_ids) - set(mapping))
    if missing_samples:
        raise ValueError("Sample(s) missing from mapping CSV: " + ", ".join(missing_samples))

    extra_samples = sorted(set(mapping) - set(sample_ids))
    if extra_samples:
        print("Warning: mapping rows not present in input: " + ", ".join(extra_samples))

    return mapping


def main():
    parser = argparse.ArgumentParser(
        description="Prepare clone-consensus input from nf-viral-integration final_results"
    )
    parser.add_argument("final_results")
    parser.add_argument("output_csv")
    parser.add_argument("--sample-map", default="-")
    args = parser.parse_args()

    final_results = Path(args.final_results)
    output_csv = Path(args.output_csv)
    sample_map = None if args.sample_map == "-" else Path(args.sample_map)

    if not final_results.is_dir():
        raise SystemExit(f"Missing final-results folder: {final_results}")
    if sample_map is not None and not sample_map.is_file():
        raise SystemExit(f"Missing sample mapping CSV: {sample_map}")

    files = sorted(final_results.glob("*/*.annotated.csv"))
    if not files:
        raise SystemExit("No *.annotated.csv files found in sample subfolders")

    sample_files = {}
    for file in files:
        sample_id = file.name[:-len(".annotated.csv")]
        if not sample_id or BAD_ID.search(sample_id):
            raise ValueError(f"Invalid sample_id from filename: {file.name}")
        if sample_id in sample_files:
            raise ValueError(f"Duplicate sample_id found: {sample_id}")
        sample_files[sample_id] = file

    mapping = load_mapping(sample_map, sample_files)
    output_rows = []
    skipped_nonflanked = 0
    skipped_incomplete = 0

    for sample_id, file in sample_files.items():
        with open(file, newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                raise ValueError(f"No header in {file}")

            required = ["READ", "STRAND", "HIV_SEQ", "INTEGRATION_SITE"]
            missing = [x for x in required if x not in reader.fieldnames]
            if missing:
                raise ValueError(f"Missing column(s) in {file}: {', '.join(missing)}")

            chromosome_col = next(
                (x for x in ["chromosome", "CHROMOSOME", "chrom"] if x in reader.fieldnames),
                None,
            )
            if chromosome_col is None:
                raise ValueError(f"Missing chromosome column in {file}")

            for row in reader:
                chromosome = clean(row[chromosome_col])
                integration_site = clean(row["INTEGRATION_SITE"])
                hiv_seq = clean(row["HIV_SEQ"])

                # Standard nf-viral-integration output provides a host integration
                # coordinate only for flanked reads.
                if chromosome.upper() == "HIV":
                    skipped_nonflanked += 1
                    continue
                if not chromosome or not integration_site or not hiv_seq:
                    skipped_incomplete += 1
                    continue

                output_rows.append({
                    "participant_id": mapping[sample_id],
                    "sample_id": sample_id,
                    "read": clean(row["READ"]),
                    "clone_id": f"{chromosome}_{integration_site}",
                    "STRAND": clean(row["STRAND"]).lower(),
                    "HIV_SEQ": hiv_seq,
                })

    if not output_rows:
        raise ValueError("No flanked HIV fragments with clone assignments were found")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REQUIRED_OUTPUT)
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Prepared {len(output_rows)} flanked fragments from {len(sample_files)} samples")
    print(f"Skipped non-flanked fragments: {skipped_nonflanked}")
    if skipped_incomplete:
        print(f"Skipped incomplete flanked rows: {skipped_incomplete}")
    print(f"Clone-consensus input: {output_csv}")


if __name__ == "__main__":
    main()
