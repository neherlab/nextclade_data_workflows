import numpy as np
import pandas as pd
import typer
from Bio import SeqIO
from joblib import Parallel, cpu_count, delayed
from tqdm import tqdm  # Optional: for progress bar

# --- Global variables for worker processes ---
# These will be initialized in main() and inherited by worker processes.
_ref_vec_worker = None
_reference_length_worker = None


# --- Helper functions for parallel processing (must be top-level) ---
def char_to_int_worker(char):
    if char == "A":
        return 0
    elif char == "C":
        return 1
    elif char == "G":
        return 2
    return 3  # T or other


def process_ranges_worker(vec, ranges_str, index_to_set_at_one):
    if pd.isna(ranges_str) or not ranges_str.strip():
        return vec
    for i in ranges_str.split(","):
        item = i.strip()
        if not item:
            continue
        split = item.split("-")
        start_pos, end_pos = 0, 0
        try:
            if len(split) == 1:
                start_pos = end_pos = int(split[0])
            elif len(split) == 2:
                start_pos = int(split[0])
                end_pos = int(split[1])
            else:
                continue  # Malformed range
        except ValueError:
            # print(f"Warning: Malformed range item '{item}'")
            continue

        # Assuming positions are 1-based, convert to 0-based index for vec
        # Loop from start_pos-1 up to end_pos-1
        for j in range(start_pos - 1, end_pos):
            if 0 <= j < _reference_length_worker:  # Check bounds
                base_idx = 6 * j
                if base_idx + 5 < len(vec):  # Ensure full segment exists
                    vec[base_idx : base_idx + 5] = 0  # Clear ACGTD flags
                    vec[base_idx + index_to_set_at_one] = 1
    return vec


def process_missing_worker(vec, missing_str):
    return process_ranges_worker(vec, missing_str, 5)  # 5 for 'missing' flag


def process_deletions_worker(vec, deletions_str):
    return process_ranges_worker(vec, deletions_str, 4)  # 4 for 'deletion' flag


def process_substitutions_worker(vec, subs_str):
    if pd.isna(subs_str) or not subs_str.strip():
        return vec
    for sub in subs_str.split(","):
        sub_item = sub.strip()
        if not sub_item:
            continue
        try:
            # Assuming format like C123A (RefChar Pos AltChar)
            # Original position is 1-based
            pos = int(sub_item[1:-1]) - 1  # Convert to 0-based index
            alt_char = sub_item[-1]
            char_idx = char_to_int_worker(alt_char)

            if 0 <= pos < _reference_length_worker:  # Check bounds
                base_idx = 6 * pos
                if base_idx + 5 < len(vec):  # Ensure full segment exists
                    vec[base_idx : base_idx + 5] = 0  # Clear ACGTD flags
                    vec[base_idx + char_idx] = 1
        except (ValueError, IndexError):
            # print(f"Warning: Could not parse substitution '{sub_item}'")
            continue
    return vec


def process_row_worker(row_missing, row_deletions, row_substitutions):
    # _ref_vec_worker is a global numpy array in the worker process
    v = _ref_vec_worker.copy()
    process_missing_worker(v, row_missing)
    process_deletions_worker(v, row_deletions)
    process_substitutions_worker(v, row_substitutions)
    return v


def process_group_worker(group_tuple):
    # group_tuple is (name, group_df)
    # _reference_length_worker is a global int in the worker process
    name, group_df = group_tuple
    # Initialize accumulator for this group
    acc = np.zeros(6 * _reference_length_worker, dtype=np.int32)
    for _, row in group_df.iterrows():
        # Pass relevant Series values to process_row_worker
        acc += process_row_worker(
            row["missing"], row["deletions"], row["substitutions"]
        )
    return name, acc


# --- append_tails function (used sequentially before parallel part) ---
def append_tails_corrected(
    alignment_start,
    alignment_end,
    current_missing_str,
    non_acgtns_str,
    total_ref_length,
):
    missing_parts = []
    if pd.notna(current_missing_str) and str(current_missing_str).strip() != "":
        missing_parts.append(str(current_missing_str).strip())

    start, end = pd.NA, pd.NA
    if pd.notna(alignment_start):
        try:
            start = int(alignment_start)
        except ValueError:
            pass
    if pd.notna(alignment_end):
        try:
            end = int(alignment_end)
        except ValueError:
            pass

    if pd.notna(start) and start > 1:
        missing_parts.append(f"1-{start - 1}")

    if pd.notna(end) and end < total_ref_length:
        missing_parts.append(f"{end + 1}-{total_ref_length}")

    if pd.notna(non_acgtns_str) and isinstance(non_acgtns_str, str):
        ambs = non_acgtns_str.split(",")
        for amb_entry in ambs:
            entry = amb_entry.strip()
            if not entry:
                continue
            try:
                parts = entry.split(":")
                # Expecting format like "TYPE:POSITION" or "REF:POSITION:ALT"
                if len(parts) > 1 and parts[1].isdigit():
                    missing_parts.append(parts[1])
                # else:
                #     print(f"Warning: Could not extract numeric position from ambiguous entry '{entry}'")
            except Exception:
                # print(f"Warning: Error parsing ambiguous entry '{entry}'")
                continue

    final_missing_str = ",".join(filter(None, missing_parts))
    return final_missing_str.replace(",,", ",").strip(",")


def main(
    ref: str = "reference.fasta",
    meta: str = "meta.tsv",
    out: str = "pango_matrix.npz",
):
    global _ref_vec_worker, _reference_length_worker

    df = pd.read_csv(meta, sep="\t")
    reference_sequence = str(SeqIO.read(ref, "fasta").seq)

    # Initialize global variables for workers
    _reference_length_worker = len(reference_sequence)

    _ref_vec_worker = np.zeros(6 * _reference_length_worker, dtype=np.int32)
    for i, c in enumerate(reference_sequence):
        _ref_vec_worker[i * 6 + char_to_int_worker(c)] = 1

    # Sequential step: Modify 'missing' column
    # Ensure relevant columns exist, use .get() for safety if columns might be missing.
    df["missing"] = df.apply(
        lambda row: append_tails_corrected(
            row.get("alignmentStart"),
            row.get("alignmentEnd"),
            row.get("missing"),
            row.get("nonACGTNs"),
            _reference_length_worker,
        ),
        axis=1,
    )

    # Prepare groups for parallel processing
    # df.groupby() returns a generator. Convert to list for joblib.
    grouped_data = list(df.groupby("pango_designated"))

    num_cores = cpu_count()  # Use all available cores
    print(f"Processing {len(grouped_data)} groups using {num_cores} cores...")

    # Parallel execution using joblib
    # The `delayed` function wraps the function and its arguments for lazy evaluation.
    # `Parallel` executes these delayed calls in parallel.
    # tqdm can be wrapped around `grouped_data` for progress display.
    results_list = Parallel(n_jobs=num_cores)(
        delayed(process_group_worker)(item)
        for item in tqdm(grouped_data, desc="Processing groups")
    )

    # Convert list of (name, acc_vector) tuples back to a dictionary
    result_dict = dict(results_list)

    np.savez_compressed(out, **result_dict)
    print(f"Output saved to {out}")


if __name__ == "__main__":
    typer.run(main)
