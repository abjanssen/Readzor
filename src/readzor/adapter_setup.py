"""Adapter loading and per-mate adapter selection."""

import os
import re

from readzor.constants import DEFAULT_ADAPTERS


def load_adapters_from_arg(value):
    """
    Load adapter sequences given on the command line, either as a FASTA
    file path or as comma-separated sequences.

    If ``value`` is an existing file, it is parsed with
    load_adapters_from_fasta(). Otherwise it is split on commas; each
    sequence is stripped, uppercased and must consist of A, C, G, T and N
    only. Such sequences are named "custom_1", "custom_2", ...

    Args:
        value (str): FASTA file path or comma-separated sequences.

    Returns:
        list[tuple[str, str]]: (name, sequence) tuples.

    Raises:
        ValueError: If the FASTA file holds no sequences, no sequences are
            given, or an entry is neither an existing file nor a valid
            sequence.
    """
    if os.path.isfile(value):
        adapters = load_adapters_from_fasta(value)
        if not adapters:
            raise ValueError(f"No sequences in file '{value}' detected.")
        return adapters
    seqs = [s.strip().upper() for s in value.split(",") if s.strip()]
    if not seqs:
        raise ValueError("No adapter sequences detected.")
    for s in seqs:
        if not re.fullmatch(r"[ACGTN]+", s):
            raise ValueError(f"'{s}' is neither an existing FASTA file nor a valid sequence.")
    return [(f"custom_{i}", s) for i, s in enumerate(seqs, 1)]


def load_adapters_from_fasta(fasta_file):
    """
    Parse adapter sequences from a FASTA file.

    Each entry is a header line starting with ">" followed by one or more
    sequence lines. Sequence lines are uppercased and concatenated, so a
    sequence may be split over several lines. Trailing whitespace is
    stripped from every line.

    Note: parsing stops at the first empty (or whitespace-only) line,
    including one at the very top of the file, so the file should not
    contain blank lines.

    Args:
        fasta_file (str): Path to the FASTA file containing adapter sequences.

    Returns:
        list[tuple[str, str]]: A list of (name, sequence) tuples, where:
            - name (str): The header line without the leading ">" and
              surrounding whitespace.
            - sequence (str): The full, uppercased nucleotide sequence.

    Raises:
        ValueError: If an entry does not start with a ">" header line, or a
            header is not followed by any sequence lines.
    """
    result = []
    with open(fasta_file) as fastafile:
        header = fastafile.readline().rstrip()
        while header:
            if not header.strip():
                header = fastafile.readline().rstrip()
                continue
            if not header.startswith(">"):
                raise ValueError("Malformed adapter Fasta file detected.")
            output_sequence = []
            sequence = fastafile.readline().rstrip()
            while sequence and not sequence.startswith(">"):
                output_sequence.append(sequence.upper())
                sequence = fastafile.readline().rstrip()
            if len(output_sequence) == 0:
                raise ValueError(f"No sequence detected for header {header}")
            name = header.lstrip(">").strip()
            joined_sequence = "".join(output_sequence)
            result.append((name, joined_sequence))
            header = sequence
        return result


def select_adapters_for_read(adapter_sequences, adapter_group, read_direction,adapter_per_read=None):
    """
    Select the adapter sequences to search for in one mate of a read pair.

    If mate-specific adapters were given (--adapter-read-1/--adapter-read-2,
    passed as ``adapter_per_read``), the list for ``read_direction`` is
    returned as is; unpaired reads get all of ``adapter_sequences``.

    Otherwise, for paired reads (read_direction "read_1" or "read_2"), the
    other mate's directional adapter (a DEFAULT_ADAPTERS entry named
    "Read 1" or "Read 2") of each selected adapter group is dropped, e.g.
    the TruSeq "Read 2" adapter is not searched for in R1. A dropped
    sequence is kept anyway if another entry of a selected group uses the
    same sequence. Other sequences (e.g. from --adapter-add) are kept
    unless they are identical to a dropped one.

    Args:
        adapter_sequences (list[bytes]): All adapter sequences of the run.
        adapter_group (list[str] | None): Selected built-in adapter groups,
            or None (as with --adapter-excl) to skip the selection.
        read_direction (str | None): "read_1", "read_2", or None for
            unpaired reads.
        adapter_per_read (dict[str, list[bytes]] | None): Mate-specific
            adapters keyed "read_1"/"read_2", or None. Defaults to None.

    Returns:
        list[bytes]: The adapter sequences to search for. Without
            ``adapter_per_read``, returned unchanged if read_direction is
            not "read_1"/"read_2" or adapter_group is empty or None.
    """
    if adapter_per_read:
        if read_direction in adapter_per_read:
            return adapter_per_read[read_direction]
        return adapter_sequences
    
    if read_direction not in ("read_1", "read_2") or not adapter_group:
        return adapter_sequences
    other_name = "Read 2" if read_direction == "read_1" else "Read 1"
    drop, keep = set(), set()
    for group, entries in DEFAULT_ADAPTERS:
        if group not in adapter_group:
            continue
        for name, seq in entries:
            data = seq.encode("utf-8")
            if name == other_name:
                drop.add(data)
            else:
                keep.add(data)
    return [a for a in adapter_sequences if a not in drop or a in keep]
