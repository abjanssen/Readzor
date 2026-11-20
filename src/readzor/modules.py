"""Trimming and filtering modules.

Each module takes per-chunk arrays and returns per-read (left, right)
boundaries of the bases it would keep."""

import re

import numpy as np

from readzor.constants import COMP_TABLE
from readzor.adapter_setup import select_adapters_for_read
from readzor.arrays import average_quality_batch

def kmer_complexity_scan(sequence_arr, chunk_padding_bool, padding_mask_bool, kmer, low_complex_cutoff, allow_n):
    """
    Low-complexity filter: flag reads with too few distinct k-mers for removal.

    For each k-mer length, the number of distinct k-mers in a read is divided
    by the maximum possible number, min(number of k-mer windows,
    alphabet_size ** k). A read fails if this ratio is below
    `low_complex_cutoff` percent for any of the k-mer lengths. Failing reads
    get an empty window; passing reads are not trimmed.

    Bases are encoded with 3 bits each, which limits k to 21 (63 bits). The
    alphabet size is 5 (A, C, G, T, N) when `allow_n` is True, and 4
    otherwise; N is then not a separate letter (it is encoded like A), as
    reads with N are expected to have been removed by --n-filter. In padded
    chunks, k-mer windows overlapping padding are excluded from both the
    distinct count and the maximum.

    Args:
        sequence_arr (numpy.ndarray): A 2D array of ASCII sequence codes
            of shape (n_reads, length).
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.
        padding_mask_bool (numpy.ndarray | None): (n_reads, length) boolean
            mask marking right-padding positions, or None if the chunk
            isn't padded.
        kmer (int | str | iterable): The k-mer length(s) to evaluate. A
            string may contain comma-separated values, e.g. "4,6".
        low_complex_cutoff (float): Minimum percentage of distinct k-mers,
            relative to the maximum possible, required to pass.
        allow_n (bool): If True, N is counted as a fifth letter of the
            alphabet.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,). Left cutoffs are all 0 (int8); right
            cutoffs (int32) are the chunk's row width for passing reads and
            0 for low-complexity reads.

    Raises:
        ValueError: If any k-mer length is greater than the chunk's row
            width, or greater than 21.
    """
    n_reads, length = sequence_arr.shape
    if isinstance(kmer, str):
        kmer_list = [int(k.strip()) for k in kmer.split(',')]
    elif hasattr(kmer, '__iter__'):
        kmer_list = [int(k) for k in kmer]
    else:
        kmer_list = [int(kmer)]

    for k in kmer_list:
        if k > length:
            raise ValueError(f"k-mer length {k} is greater than sequence length {length}")
        if k > 21:
            raise ValueError(
                f"Chosen k-mer length {k} is too big, maximum safe k-mer is 21. "
            )

    global_passed = np.ones(n_reads, dtype=bool)
    mapping = np.zeros(256, dtype=np.int8)
    if allow_n:
        mapping[ord('A')] = 0
        mapping[ord('C')] = 1
        mapping[ord('G')] = 2
        mapping[ord('T')] = 3
        mapping[ord('N')] = 4
        mapping[ord('\x00')] = 5
        bits_per_base = 3
    else:
        mapping[ord('A')] = 0
        mapping[ord('C')] = 1
        mapping[ord('G')] = 2
        mapping[ord('T')] = 3
        mapping[ord('\x00')] = 4
        bits_per_base = 3
    int_matrix_full = mapping[sequence_arr]

    for k in kmer_list:
        if not np.any(global_passed):
            break
        if k > length:
            raise ValueError(f"k-mer length {k} is greater than sequence length {length}")
            
        max_kmers = length - k + 1
        alphabet_size = 5 if allow_n else 4
        max_possible_unique = min(max_kmers, alphabet_size**k)
        
        n_bits = bits_per_base * k
        kmer_dtype = np.int16 if n_bits <= 15 else (np.int32 if n_bits <= 31 else np.int64)
        kmer_ints = np.zeros((n_reads, max_kmers), dtype=kmer_dtype)
        
        if not chunk_padding_bool:
            for i in range(k):
                kmer_ints = (kmer_ints << bits_per_base) | int_matrix_full[:, i:i+max_kmers]

            sorted_kmers = np.sort(kmer_ints, axis=1)
            is_new = np.empty_like(sorted_kmers, dtype=bool)
            is_new[:, 0] = True
            np.not_equal(sorted_kmers[:, 1:], sorted_kmers[:, :-1], out=is_new[:, 1:])
            unique_counts = is_new.sum(axis=1)

            ratio = unique_counts / max_possible_unique
            global_passed &= (ratio >= (low_complex_cutoff / 100))
        else:
            window_has_pad = np.zeros((n_reads, max_kmers), dtype=bool)
            for i in range(k):
                kmer_ints = (kmer_ints << bits_per_base) | int_matrix_full[:, i:i+max_kmers]
                window_has_pad |= padding_mask_bool[:, i:i+max_kmers]
            kmer_ints = np.where(window_has_pad, -1, kmer_ints)

            sorted_kmers = np.sort(kmer_ints, axis=1)
            is_new = np.empty_like(sorted_kmers, dtype=bool)
            is_new[:, 0] = True
            np.not_equal(sorted_kmers[:, 1:], sorted_kmers[:, :-1], out=is_new[:, 1:])
            unique_counts = is_new.sum(axis=1)

            valid_kmer_counts = (~window_has_pad).sum(axis=1)
            any_pad_in_row = window_has_pad.any(axis=1)
            unique_counts = unique_counts - any_pad_in_row.astype(np.int32)

            denom = np.minimum(valid_kmer_counts, alphabet_size**k)

            ratio = np.where(denom > 0, unique_counts / denom, 0.0)
            global_passed &= (ratio >= (low_complex_cutoff / 100))

    second_array = np.where(global_passed, length, 0).astype(np.int32)
    return np.zeros(n_reads, dtype=np.int8), second_array

def homopolymer_nucleotide_trimming(sequence_arr, padding_mask_bool, row_tilde_count, chunk_padding_bool, poly_length_start, poly_length_end, poly_bases_start, poly_bases_end):
    """
    Homopolymer trimming: per-read boundaries that remove homopolymer runs
    from the start and/or end of each read.

    For each base to check at an end, the run of that base touching the end
    is measured; if it is at least the minimum length for that end, the whole
    run is trimmed. Each listed base is checked on its own (a mixed run such
    as "GGAGG" is not treated as one run), and the largest trim per end wins.

    Each end is configured separately (the --poly-bases-both and
    --poly-length-both fallbacks are already resolved in parse_args()). An
    end is only checked if it has at least one base and a minimum length
    greater than 0; a length of None counts as 0. If neither end is
    checked, no trimming is done. A read consisting entirely of a checked
    base ends up with an empty window.

    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes.
        padding_mask_bool (numpy.ndarray | None): (n_reads, read_length)
            boolean mask marking right-padding positions, or None if the
            chunk isn't padded.
        row_tilde_count (numpy.ndarray | int): (n_reads,) count of padding
            bytes per read, used when chunk_padding_bool is True.
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.
        poly_length_start (int | None): Minimum run length to trigger
            trimming at the start of the read (0 or None = start not checked).
        poly_length_end (int | None): Minimum run length to trigger
            trimming at the end of the read (0 or None = end not checked).
        poly_bases_start (str | None): Comma-separated bases to check at
            the start (case-insensitive), or None for none.
        poly_bases_end (str | None): Comma-separated bases to check at
            the end (case-insensitive), or None for none.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,), giving the left and right trim
            boundaries per read.

    Raises:
        ValueError: If any comma-separated base entry is not a single character.
    """
    n_reads, length = sequence_arr.shape
    
    def parse_bases(value):
        if not value:
            return []
        bases = [b.strip().upper() for b in value.split(",") if b.strip()]
        if any(len(b) != 1 for b in bases):
            raise ValueError(f"Invalid base entry in '{value}': All bases must be single characters.")
        return bases
    
    start_bases = parse_bases(poly_bases_start)
    end_bases = parse_bases(poly_bases_end)
    if poly_length_start is None:
        poly_length_start = 0
    if poly_length_end is None:
        poly_length_end = 0
        
    if poly_length_start <= 0:
        start_bases = []
    if poly_length_end <= 0:
        end_bases = []
    if not start_bases and not end_bases:
        return np.zeros(n_reads, dtype=np.int8), np.full(n_reads, length, dtype=np.int32)
    
    right_cutoffs = np.full(n_reads, length, dtype=np.int32)
    left_cutoffs = np.zeros(n_reads, dtype=np.int8)
    for base in start_bases:
        base_code = ord(base)
        non_base_mask = sequence_arr != base_code
        first_non_pos = non_base_mask.argmax(axis=1)
        first_non_pos = np.where(non_base_mask.any(axis=1), first_non_pos, length)
        trim_amount = np.where(first_non_pos >= poly_length_start, first_non_pos, 0)
        left_cutoffs = np.maximum(left_cutoffs, trim_amount)

    if end_bases:
        rev_seq = np.ascontiguousarray(sequence_arr[:, ::-1])
        if not chunk_padding_bool:
            for base in end_bases:
                base_code = ord(base)
                non_base_mask = rev_seq != base_code
                first_non_pos = non_base_mask.argmax(axis=1)
                first_non_pos = np.where(non_base_mask.any(axis=1), first_non_pos, length)
                trim_amount = np.where(first_non_pos >= poly_length_end, first_non_pos, 0)
                base_right_cutoffs = length - trim_amount
                right_cutoffs = np.minimum(right_cutoffs, base_right_cutoffs)
        else:
            real_length = length - row_tilde_count
            for base in end_bases:
                non_base_mask = (rev_seq != ord(base)) & (rev_seq != ord('\x00'))
                first_non_pos = non_base_mask.argmax(axis=1)
                first_non_pos = np.where(non_base_mask.any(axis=1), first_non_pos, length)
                run_length = first_non_pos - row_tilde_count
                trim_amount = np.where(run_length >= poly_length_end, run_length, 0)
                base_right_cutoffs = real_length - trim_amount
                right_cutoffs = np.minimum(right_cutoffs, base_right_cutoffs)

    return left_cutoffs, right_cutoffs

def adapter_trimming(sequence_arr, chunk_padding_bool, row_tilde_count, adapter_sequences, mismatches, adapter_seed, adapter_per_read, read_direction=None, adapter_group=None):
    """
    Adapter trimming: cut each read at the earliest occurrence of any
    adapter sequence.

    For each read, the earliest position where any of the adapters occurs is
    found, and everything from that position onwards is removed (3'-end
    trimming). An adapter that runs off the 3' end of the read is also
    detected, as long as at least `adapter_seed` of its bases (or the whole
    adapter, if it is shorter) lie within the read. For paired reads,
    overlap trimming can catch shorter adapter remnants.

    With mismatches == 0 and an effective seed length (the smaller of
    `adapter_seed` and the shortest adapter) of at least 4, matching is
    exact: the adapter's first `adapter_seed` bases are located with a
    substring search, and each hit is then checked against the rest of the
    adapter, or against the part of it that fits in the read. Otherwise, a
    vectorized Hamming-distance search is used, which allows substitutions
    (no insertions or deletions). The allowed number scales with the
    overlap: L adapter bases within the read may contain
    floor(L * mismatches / adapter_length) mismatches.

    For paired reads, the adapters searched for can be narrowed down per
    mate first (see select_adapters_for_read()).

    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes (int8).
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.
        row_tilde_count (numpy.ndarray | int): (n_reads,) count of padding
            bytes per read, used to keep matches within each read's real
            length when chunk_padding_bool is True.
        adapter_sequences (list[bytes]): Adapter sequences to search for.
        mismatches (int): Number of substitutions allowed over the full
            adapter length. 0 means exact matching.
        adapter_seed (int): Minimum number of adapter bases that must lie
            within the read. With exact matching, also the length of the
            adapter prefix used for the substring search.
        adapter_per_read (dict[str, list[bytes]] | None): Mate-specific
            adapters keyed "read_1"/"read_2" (from --adapter-read-1/
            --adapter-read-2), or None.
        read_direction (str | None): "read_1" or "read_2" for one mate of a
            read pair, or None for unpaired reads. Defaults to None.
        adapter_group (list[str] | None): Selected built-in adapter groups,
            used together with read_direction to select the adapters per
            mate, or None to search for all of adapter_sequences.
            Defaults to None.
    Returns:
        tuple[numpy.ndarray, numpy.ndarray, dict]: (left_cutoffs,
            right_cutoffs, hits), the first two of shape (n_reads,). Left
            cutoffs are always 0 (int8); right cutoffs (int32) are the
            adapter start position, or the read's real length if no adapter
            was found. ``hits`` holds "adapters" (the sequences searched for,
            after per-mate selection) and "hit_index" ((n_reads,) int32:
            index into "adapters" of the adapter that set the cut, or -1),
            used by --stats.
    """
    n_reads, length = sequence_arr.shape
    adapter_sequences = select_adapters_for_read(adapter_sequences, adapter_group, read_direction, adapter_per_read)
    if not adapter_sequences:
        if chunk_padding_bool:
            right_cutoffs = (length - row_tilde_count).astype(np.int32)
        else:
            right_cutoffs = np.full(n_reads, length, dtype=np.int32)
        return np.zeros(n_reads, dtype=np.int8), right_cutoffs, {"adapters": adapter_sequences, "hit_index": np.full(n_reads, -1, dtype=np.int32)}
    hit_index = np.full(n_reads, -1, dtype=np.int32)
    adapter_len_min = min(len(seq) for seq in adapter_sequences)
    if mismatches == 0 and min(adapter_seed, adapter_len_min) >= 4:
        all_bytes = sequence_arr.tobytes()
        if not chunk_padding_bool:
            real_lengths = np.full(n_reads, length, dtype=np.int32)
        else:
            real_lengths = length - row_tilde_count
        right_cutoffs = np.array(real_lengths, dtype=np.int32)
        for adapter_idx, adapter_bytes in enumerate(adapter_sequences):
            adapter_len = len(adapter_bytes)
            seed_len = min(adapter_seed, adapter_len)
            seed = adapter_bytes[:seed_len]
            start = 0
            while True:
                pos = all_bytes.find(seed, start)
                if pos == -1:
                    break
                i, pos_in_row = divmod(pos, length)
                avail = real_lengths[i] - pos_in_row
                if avail >= adapter_len:
                    if all_bytes.startswith(adapter_bytes, pos):
                        if pos_in_row < right_cutoffs[i]:
                            right_cutoffs[i] = pos_in_row
                            hit_index[i] = adapter_idx
                        start = (i + 1) * length
                    else:
                        start = pos + 1
                    continue
                if avail == seed_len:
                    if pos_in_row < right_cutoffs[i]:
                        right_cutoffs[i] = pos_in_row
                        hit_index[i] = adapter_idx
                    start = (i + 1) * length
                    continue
                if avail > seed_len:
                    if all_bytes[pos:pos + avail] == adapter_bytes[:avail]:
                        if pos_in_row < right_cutoffs[i]:
                            right_cutoffs[i] = pos_in_row
                            hit_index[i] = adapter_idx
                        start = (i + 1) * length
                    else:
                        start = pos + 1
                    continue
                start = pos + 1
    else:
        if not chunk_padding_bool:
            right_cutoffs = np.full(n_reads, length, dtype=np.int32)
            pos = np.arange(length, dtype=np.int32)
            for adapter_idx, adapter_bytes in enumerate(adapter_sequences):
                adapter_len = len(adapter_bytes)
                seed_len = min(adapter_seed, adapter_len)
                n_steps = min(adapter_len, length)
                mm_dtype = np.int8 if (n_steps <= 127 and mismatches <= 127) else np.int16

                n_cols = length - seed_len + 1
                if n_cols <= 0:
                    continue
                overlap = np.minimum(adapter_len, length - pos[:n_cols])
                allowed = np.where(overlap >= seed_len, (overlap * mismatches) // adapter_len, -1).astype(mm_dtype)
                check_step = 14 + 2 * mismatches
                n_open = min(length - check_step, n_cols)
                do_check = check_step < n_steps and n_open > 0

                mismatch_matrix = np.zeros((n_reads, n_cols), dtype=mm_dtype)

                first = check_step if do_check else n_steps
                for j in range(first):
                    w = min(length - j, n_cols)
                    mismatch_matrix[:, :w] += sequence_arr[:, j:j + w] != adapter_bytes[j]

                if do_check:
                    alive = (mismatch_matrix[:, :n_open] <= allowed[:n_open]).any(axis=1)
                    valid_mask = mismatch_matrix <= allowed
                    valid_mask[~alive, :n_open] = False
                    idx = np.flatnonzero(alive)
                    if idx.size:
                        sub = mismatch_matrix[idx]
                        seq_sub = sequence_arr[idx]
                        for j in range(first, n_steps):
                            w = min(length - j, n_cols)
                            sub[:, :w] += seq_sub[:, j:j + w] != adapter_bytes[j]
                        valid_mask[idx] = sub <= allowed
                else:
                    valid_mask = mismatch_matrix <= allowed

                has_match = valid_mask.any(axis=1)
                if has_match.any():
                    first_match_col = valid_mask.argmax(axis=1)
                    improved = has_match & (first_match_col < right_cutoffs)
                    right_cutoffs[improved] = first_match_col[improved].astype(np.int32)
                    hit_index[improved] = adapter_idx
        else:
            max_adapter_len = max(len(a) for a in adapter_sequences)
            idx_dtype = np.int16 if max(length + max_adapter_len, max_adapter_len * mismatches) < 32000 else np.int32

            real_lengths = (length - row_tilde_count).astype(idx_dtype)
            right_cutoffs = real_lengths.astype(np.int32)
            pos = np.arange(length, dtype=idx_dtype)
            rl = real_lengths[:, None]
            for adapter_idx, adapter_bytes in enumerate(adapter_sequences):
                adapter_len = len(adapter_bytes)
                seed_len = min(adapter_seed, adapter_len)
                n_steps = min(adapter_len, length)
                mm_dtype = np.int8 if (n_steps <= 127 and mismatches <= 127) else np.int16
                n_cols = length - seed_len + 1
                if n_cols <= 0:
                    continue
                pos_c = pos[:n_cols]
                overlap = np.minimum(adapter_len, rl - pos_c)
                seed_ok = overlap >= seed_len
                allowed = ((np.maximum(overlap, 0) * mismatches) // adapter_len).astype(mm_dtype)

                check_step = 14 + 2 * mismatches
                n_open = min(length - check_step, n_cols)
                do_check = check_step < n_steps and n_open > 0

                mismatch_matrix = np.zeros((n_reads, n_cols), dtype=mm_dtype)

                first = check_step if do_check else n_steps
                for j in range(first):
                    w = min(length - j, n_cols)
                    mismatch_matrix[:, :w] += sequence_arr[:, j:j + w] != adapter_bytes[j]

                if do_check:
                    pad_partial = np.clip(np.minimum(pos_c + check_step, length) - rl, 0, n_steps).astype(mm_dtype)
                    corrected = mismatch_matrix - pad_partial
                    valid_mask = seed_ok & (corrected <= allowed)
                    alive = valid_mask[:, :n_open].any(axis=1)
                    valid_mask[~alive, :n_open] = False
                    idx = np.flatnonzero(alive)
                    if idx.size:
                        sub = mismatch_matrix[idx]
                        seq_sub = sequence_arr[idx]
                        for j in range(first, n_steps):
                            w = min(length - j, n_cols)
                            sub[:, :w] += seq_sub[:, j:j + w] != adapter_bytes[j]
                        sub -= np.clip(np.minimum(pos_c + adapter_len, length) - rl[idx], 0, n_steps).astype(mm_dtype)
                        valid_mask[idx] = seed_ok[idx] & (sub <= allowed[idx])
                else:
                    mismatch_matrix -= np.clip(np.minimum(pos_c + adapter_len, length) - rl, 0, n_steps).astype(mm_dtype)
                    valid_mask = seed_ok & (mismatch_matrix <= allowed)

                has_match = valid_mask.any(axis=1)
                if has_match.any():
                    first_match_col = valid_mask.argmax(axis=1)
                    improved = has_match & (first_match_col < right_cutoffs)
                    right_cutoffs[improved] = first_match_col[improved].astype(np.int32)
                    hit_index[improved] = adapter_idx
    return np.zeros(n_reads, dtype=np.int8), right_cutoffs, {"adapters": adapter_sequences, "hit_index": hit_index}

def adapter_trimming_overlap(seq_arr_1, seq_arr_2, chunk_padding_bool_1, row_tilde_count_1, chunk_padding_bool_2, row_tilde_count_2, min_overlap, max_mismatch, max_mismatch_frac):
    """
    Overlap-based adapter trimming for read pairs, independent of the
    adapter sequence.

    When the DNA insert is shorter than the reads, both mates read through
    into adapter. For a pair with insert size I, the first I bases of R1
    then equal the reverse complement of the first I bases of R2. For every
    candidate I from `min_overlap` upwards (only sizes that both mates cover
    and that at least one mate reads past), these two I-base regions are
    compared. Any position where either base is N counts as a mismatch,
    including N vs N. An insert size is accepted if it has at most
    `max_mismatch` mismatches and a mismatch fraction of at most
    `max_mismatch_frac` percent of the I overlap positions. If several
    insert sizes are accepted, the largest one wins. Both mates are then
    cut at I.

    As a speed-up, each candidate is first compared on its first
    2 * max_mismatch + 6 positions only; candidates that already exceed
    `max_mismatch` there are skipped.

    Args:
        seq_arr_1 (numpy.ndarray): (n_reads, length_1) int8 ASCII array of R1.
        seq_arr_2 (numpy.ndarray): (n_reads, length_2) int8 ASCII array of R2,
            row-aligned with seq_arr_1 (row i of both is one pair).
        chunk_padding_bool_1 (bool): Whether seq_arr_1 contains padded
            (unequal-length) reads.
        row_tilde_count_1 (numpy.ndarray | int): (n_reads,) count of padding
            bytes per R1 read, used when chunk_padding_bool_1 is True.
        chunk_padding_bool_2 (bool): As chunk_padding_bool_1, for R2.
        row_tilde_count_2 (numpy.ndarray | int): As row_tilde_count_1, for R2.
        min_overlap (int): Shortest insert size tested.
        max_mismatch (int): Maximum mismatches allowed in the overlap
            (N counts as a mismatch).
        max_mismatch_frac (float): Maximum percentage (0-100) of overlap
            positions that may mismatch (N counts as a mismatch).

    Returns:
        tuple[tuple[numpy.ndarray, numpy.ndarray], tuple[numpy.ndarray, numpy.ndarray]]:
            ((left_cutoffs_1, right_cutoffs_1), (left_cutoffs_2, right_cutoffs_2)),
            each of shape (n_reads,). Left cutoffs are always 0 (3'-end trimming
            only, dtype int8); right cutoffs (int32) are the accepted insert
            size, or the read's real length if no overlap was accepted.

    Raises:
        ValueError: If seq_arr_1 and seq_arr_2 have different numbers of rows.
    """
    prefix_length = 2 * max_mismatch + 6
    max_mismatch_frac = max_mismatch_frac

    n_base = np.int8(ord("N"))

    n_reads_1, length_1 = seq_arr_1.shape
    n_reads_2, length_2 = seq_arr_2.shape
    if n_reads_1 != n_reads_2:
        raise ValueError(f"R1 has {n_reads_1} reads, R2 has {n_reads_2}")

    rc_seq_arr_2 = np.frombuffer(seq_arr_2[:, ::-1].tobytes().translate(COMP_TABLE), dtype=np.int8).reshape(n_reads_2, length_2)    
    seq_arr_1 = np.where(seq_arr_1 == n_base, np.int8(-1), seq_arr_1)
    rc_seq_arr_2 = np.where(rc_seq_arr_2 == n_base, np.int8(-2), rc_seq_arr_2)
    
    best_insert = np.full(n_reads_1, -1, dtype=np.int32)
    max_insert_size = min(length_1, length_2)
    any_padding = chunk_padding_bool_1 or chunk_padding_bool_2

    real_lengths_1 = (length_1 - row_tilde_count_1 if chunk_padding_bool_1
                      else np.full(n_reads_1, length_1)).astype(np.int32)
    real_lengths_2 = (length_2 - row_tilde_count_2 if chunk_padding_bool_2
                      else np.full(n_reads_2, length_2)).astype(np.int32)

    if any_padding:
        shorter_mate = np.minimum(real_lengths_1, real_lengths_2)
        longer_mate = np.maximum(real_lengths_1, real_lengths_2)
        last_insert = int(np.minimum(shorter_mate, longer_mate - 1).max())
    else:
        last_insert = min(max_insert_size, max(length_1, length_2) - 1)

    for insert_size in range(min_overlap, last_insert + 1):
        offset = length_2 - insert_size
        prefix = min(prefix_length, insert_size)
        
        if any_padding:
            candidate = (insert_size <= shorter_mate) & (insert_size < longer_mate)
            candidate_rows = np.flatnonzero(candidate)
        else:
            candidate_rows = None

        if candidate_rows is None:
            prefix_arr_1 = seq_arr_1[:, :prefix]
            prefix_arr_2 = rc_seq_arr_2[:, offset:offset + prefix]
        else:
            prefix_arr_1 = seq_arr_1[candidate_rows, :prefix]
            prefix_arr_2 = rc_seq_arr_2[candidate_rows, offset:offset + prefix]
        prefix_mismatches = (prefix_arr_1 != prefix_arr_2).sum(axis=1)

        survive_bool = prefix_mismatches <= max_mismatch
        if not survive_bool.any():
            continue
        rows_survived = np.flatnonzero(survive_bool) if candidate_rows is None else candidate_rows[survive_bool]
        mismatches = prefix_mismatches[survive_bool]

        if insert_size > prefix:
            rest_arr_1 = seq_arr_1[rows_survived, prefix:insert_size]
            rest_arr_2 = rc_seq_arr_2[rows_survived, offset + prefix:]
            mismatches = mismatches + (rest_arr_1 != rest_arr_2).sum(axis=1)

        pass_bool = (mismatches <= max_mismatch) & (mismatches * 100 <= insert_size * max_mismatch_frac)
        rows_passed = rows_survived[pass_bool]
        best_insert[rows_passed] = insert_size

    found = best_insert > 0
    right_cutoffs_1 = np.where(found, best_insert, real_lengths_1).astype(np.int32)
    right_cutoffs_2 = np.where(found, best_insert, real_lengths_2).astype(np.int32)
    return ((np.zeros(n_reads_1, dtype=np.int8), right_cutoffs_1), (np.zeros(n_reads_2, dtype=np.int8), right_cutoffs_2))

def trim_ends_quality(quality_arr, chunk_padding_bool, row_tilde_count, padding_mask_bool, endqual_min_start, endqual_min_end):
    """
    Quality-dependent end trimming: per-read boundaries that trim each end
    up to the first base that meets a quality threshold.

    The left boundary is the first position (scanning from the 5' end) with
    quality >= endqual_min_start. The right boundary is one past the last
    position (scanning from the 3' end) with quality >= endqual_min_end.
    Only the ends are trimmed; low-quality bases inside the read are kept.
    Padding positions are ignored. Reads without any base meeting the
    threshold get left = row width and right = 0, i.e. an empty window.
    Vectorized across all reads in the chunk.

    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base quality scores.
        chunk_padding_bool (bool): Whether quality_arr contains padded
            (unequal-length) reads.
        row_tilde_count (numpy.ndarray | int): (n_reads,) count of padding
            bytes per read, used when chunk_padding_bool is True.
        padding_mask_bool (numpy.ndarray | None): (n_reads, read_length)
            boolean mask marking right-padding positions, or None if the
            chunk isn't padded. Padded positions are excluded from the scan.
        endqual_min_start (int): Minimum quality to keep bases from the
            start of the read. (The --endqual-min-both fallback is already
            resolved in parse_args().)
        endqual_min_end (int): Minimum quality to keep bases from the end
            of the read.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (start_cutoffs, end_cutoffs),
            each of shape (n_reads,) and dtype int32, giving the left and
            right trim boundaries per read.
    """
    n_reads, length = quality_arr.shape
    if not chunk_padding_bool:
        qual_mask = quality_arr >= endqual_min_start
        start_cutoffs = qual_mask.argmax(axis=1)
        zero_rows = start_cutoffs == 0
        if zero_rows.any():
            start_good_pos = qual_mask[:, 0] | ~zero_rows
            start_cutoffs = np.where(start_good_pos, start_cutoffs, length)

        quality_arr_rev = quality_arr[:, ::-1]
        qual_mask = quality_arr_rev >= endqual_min_end
        end_cutoffs = length - qual_mask.argmax(axis=1)
        zero_end_rows = end_cutoffs == length
        if zero_end_rows.any():
            end_good_pos = qual_mask[:, 0] | ~zero_end_rows
            end_cutoffs = np.where(end_good_pos, end_cutoffs, 0)
    else:
        qual_mask = (quality_arr >= endqual_min_start) & ~padding_mask_bool
        start_cutoffs = qual_mask.argmax(axis=1)
        zero_rows = start_cutoffs == 0
        if zero_rows.any():
            start_good_pos = qual_mask[:, 0] | ~zero_rows
            start_cutoffs = np.where(start_good_pos, start_cutoffs, length)

        quality_arr_rev = quality_arr[:, ::-1]
        pad_mask_rev = padding_mask_bool[:, ::-1]
        qual_mask = (quality_arr_rev >= endqual_min_end) & ~pad_mask_rev
        first_good_pos = qual_mask.argmax(axis=1)

        real_lengths = length - row_tilde_count
        run_length = first_good_pos - row_tilde_count
        end_cutoffs = real_lengths - run_length

        zero_end_rows = first_good_pos == 0
        if zero_end_rows.any():
            end_good_pos = qual_mask[:, 0] | ~zero_end_rows
            end_cutoffs = np.where(end_good_pos, end_cutoffs, 0)

    return start_cutoffs.astype(np.int32), end_cutoffs.astype(np.int32)

def average_quality_filter_wrapper(quality_arr, chunk_padding_bool, row_tilde_count, min_avg_qual):
    """
    Whole-read average-quality filter (--min-average-qual-pre), in the
    pipeline's (left, right) boundary format.

    The mean quality is computed over each read's full, untrimmed length.
    Reads with a mean >= min_avg_qual keep their full window
    (left 0, right = read length). Reads below the threshold get right 0,
    i.e. a zero-length window, so they are discarded downstream.

    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base quality scores.
        chunk_padding_bool (bool): Whether quality_arr contains padded
            (unequal-length) reads.
        row_tilde_count (numpy.ndarray | int): (n_reads,) count of padding
            bytes per read, used to recover each read's real length when
            chunk_padding_bool is True.
        min_avg_qual (float): Minimum average quality required to keep a read.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,). Left cutoffs are all 0 (int8); right
            cutoffs are int32.
    """
    n_reads, length = quality_arr.shape
    if not chunk_padding_bool:
        real_lengths = np.full(n_reads, length, dtype=np.int32)
    else:
        real_lengths = (length - row_tilde_count).astype(np.int32)
    avg_quals = average_quality_batch(quality_arr, lefts=0, rights=real_lengths)
    passed = avg_quals >= min_avg_qual
    right_cutoffs = np.where(passed, real_lengths, 0).astype(np.int32)
    return np.zeros(n_reads, dtype=np.int8), right_cutoffs

def sliding_window_quality(quality_arr, chunk_padding_bool, padding_mask_bool, row_tilde_count, slider_quality, slider_window, slider_step):
    """
    Sliding-window quality trimming: keep the best stretch of each read in
    which every window passes the quality threshold.

    Windows of `slider_window` bases start every `slider_step` positions
    (the last possible window is always included). A window fails if its
    mean quality is below `slider_quality`. Every base covered by at least
    one failing window is marked bad, and the longest run of good bases is
    kept, with ties broken by the higher mean quality. Because the kept
    stretch can lie in the middle of the read, both ends can be trimmed at
    once, and after a mid-read quality drop the longer side survives.

    Reads shorter than one window are not trimmed. In padded chunks,
    windows overlapping padding are not evaluated; instead the padding
    positions themselves are marked bad, and with slider_step > 1 the last
    full window of each read (ending at its real length) is also checked.

    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base quality scores.
        chunk_padding_bool (bool): Whether quality_arr contains padded
            (unequal-length) reads.
        padding_mask_bool (numpy.ndarray | None): (n_reads, read_length)
            boolean mask marking right-padding positions, or None if the
            chunk isn't padded.
        row_tilde_count (numpy.ndarray | int): (n_reads,) count of padding
            bytes per read, used to recover each read's real length when
            chunk_padding_bool is True.
        slider_quality (int): Minimum mean quality for a window to pass.
        slider_window (int): Number of bases per sliding window.
        slider_step (int): Step size between successive window start positions.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,), giving the kept [left, right) region
            per read. Reads without failing windows keep their full length;
            reads in which every base is bad get an empty region (0, 0).
    """
    n_reads, length = quality_arr.shape

    cs_dtype = np.int16 if (length * 93 < 32000 and slider_quality * slider_window < 32000) else np.int32
    pad_dtype = np.int16 if length < 32000 else np.int32

    if not chunk_padding_bool:
        if length < slider_window:
            return np.zeros(n_reads, dtype=np.int8), np.full(n_reads, length, dtype=np.int32)

        last_possible_start = length - slider_window
        window_starts = np.arange(0, last_possible_start + 1, slider_step)
        if window_starts[-1] != last_possible_start:
            window_starts = np.append(window_starts, last_possible_start)

        cumsum = np.cumsum(quality_arr, axis=1, dtype=cs_dtype)
        cumsum = np.concatenate([np.zeros((n_reads, 1), dtype=cs_dtype), cumsum], axis=1)
        window_sums = cumsum[:, window_starts + slider_window] - cumsum[:, window_starts]
        failed_mask = window_sums < (slider_quality * slider_window)

        bad_positions = np.zeros((n_reads, length), dtype=bool)

        if slider_step == 1:
            for offset in range(slider_window):
                bad_positions[:, offset:offset + len(window_starts)] |= failed_mask
        else:
            for j, start in enumerate(window_starts):
                bad_positions[:, start:start + slider_window] |= failed_mask[:, j:j + 1]

        real_lengths = np.full(n_reads, length, dtype=np.int32)
    else:
        real_lengths = length - row_tilde_count
        too_short = real_lengths < slider_window

        if length < slider_window:
            left_cutoffs = np.zeros(n_reads, dtype=np.int8)
            right_cutoffs = real_lengths.astype(np.int32)
            return left_cutoffs, right_cutoffs

        cumsum = np.cumsum(quality_arr, axis=1, dtype=cs_dtype)
        cumsum = np.concatenate([np.zeros((n_reads, 1), dtype=cs_dtype), cumsum], axis=1)
        pad_cumsum = np.cumsum(padding_mask_bool, axis=1, dtype=pad_dtype)
        pad_cumsum = np.concatenate([np.zeros((n_reads, 1), dtype=pad_dtype), pad_cumsum], axis=1)

        last_possible_start = length - slider_window
        window_starts = np.arange(0, last_possible_start + 1, slider_step)
        if window_starts[-1] != last_possible_start:
            window_starts = np.append(window_starts, last_possible_start)

        window_sums = cumsum[:, window_starts + slider_window] - cumsum[:, window_starts]
        window_pad_counts = pad_cumsum[:, window_starts + slider_window] - pad_cumsum[:, window_starts]
        window_has_pad = window_pad_counts > 0

        failed_mask = (window_sums < (slider_quality * slider_window)) & ~window_has_pad

        bad_positions = np.zeros((n_reads, length), dtype=bool)
        n_windows = len(window_starts)
        if slider_step == 1:
            for offset in range(slider_window):
                bad_positions[:, offset:offset + n_windows] |= failed_mask
        else:
            for j, start in enumerate(window_starts):
                bad_positions[:, start:start + slider_window] |= failed_mask[:, j:j + 1]

        if slider_step > 1:
            rows = np.arange(n_reads)
            has_window = real_lengths >= slider_window
            last_start = np.clip(real_lengths - slider_window, 0, None)
            last_sums = cumsum[rows, last_start + slider_window] - cumsum[rows, last_start]
            last_failed = has_window & (last_sums < slider_quality * slider_window)
            cols = np.arange(length)
            bad_positions |= (last_failed[:, None]
                              & (cols >= last_start[:, None])
                              & (cols < real_lengths[:, None]))

        bad_positions |= padding_mask_bool

    good_positions = ~bad_positions
    no_bad = ~bad_positions.any(axis=1)
    all_bad = bad_positions.all(axis=1)
    left_cutoffs = np.zeros(n_reads, dtype=np.int32)
    right_cutoffs = np.zeros(n_reads, dtype=np.int32)
    right_cutoffs[no_bad] = real_lengths[no_bad].astype(np.int32)
    needs_stretch_search = ~no_bad & ~all_bad
    if not needs_stretch_search.any():
        if chunk_padding_bool:
            left_cutoffs[too_short] = 0
            right_cutoffs[too_short] = real_lengths[too_short].astype(np.int32)
        return left_cutoffs, right_cutoffs

    padded = np.zeros((needs_stretch_search.sum(), length + 2), dtype=bool)
    padded[:, 1:-1] = good_positions[needs_stretch_search]
    diffs = np.diff(padded.view(np.int8), axis=1)
    local_rows, start_cols = np.where(diffs == 1)
    _, end_cols = np.where(diffs == -1)
    global_rows = np.where(needs_stretch_search)[0][local_rows]
    run_lengths = end_cols - start_cols
    
    run_sums = cumsum[global_rows, end_cols] - cumsum[global_rows, start_cols]
    run_means = run_sums / run_lengths
    order = np.lexsort((-run_means, -run_lengths, global_rows))
    sorted_rows = global_rows[order]
    first_in_group = np.concatenate([[0], np.flatnonzero(sorted_rows[1:] != sorted_rows[:-1]) + 1])
    best_idx = order[first_in_group]

    left_cutoffs[global_rows[best_idx]] = start_cols[best_idx]
    right_cutoffs[global_rows[best_idx]] = end_cols[best_idx]

    if chunk_padding_bool:
        left_cutoffs[too_short] = 0
        right_cutoffs[too_short] = real_lengths[too_short].astype(np.int32)

    return left_cutoffs, right_cutoffs

def cut_set_ends(sequence_arr, chunk_padding_bool, row_tilde_count, cut_start, cut_end):
    """
    Set-length end trimming: remove a fixed number of bases from each end of
    every read, regardless of sequence or quality (e.g. to hard-trim a known
    primer or a biased first/last few bases).

    `cut_start` and `cut_end` set the trim for their own end (the
    --cut-both fallback is already resolved in parse_args(), where an
    explicit 0 does override --cut-both). If the two cuts meet or cross,
    the read's window becomes empty and the read is discarded downstream.

    For equal-length chunks, the boundaries are clamped to
    [0, read_length]. For padded chunks, the boundaries are computed from
    each read's real length and are not clamped; a negative right boundary
    has the same effect (the read is discarded).

    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes, used only for its shape.
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.
        row_tilde_count (numpy.ndarray | int): (n_reads,) count of padding
            bytes per read, used to recover each read's real length when
            chunk_padding_bool is True.
        cut_start (int): Number of bases to trim from the 5' end.
        cut_end (int): Number of bases to trim from the 3' end.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,), giving the left and right trim
            boundaries per read.
    """
    n_reads, length = sequence_arr.shape
    if not chunk_padding_bool:
        cut_end_pos = length - cut_end
        cut_end_pos = max(cut_end_pos, 0)
        cut_start_pos = min(cut_start, cut_end_pos)
        cut_start_pos = max(cut_start_pos, 0)
        return np.full(n_reads, cut_start_pos, dtype=np.int32), np.full(n_reads, cut_end_pos, dtype=np.int32)
    else:
        real_length = length - row_tilde_count
        end_positions = real_length - cut_end
        cut_start_pos = np.where(cut_start > end_positions, end_positions, cut_start)
        return np.full(n_reads, cut_start_pos, dtype=np.int32), end_positions

def header_mgi_to_illumina(mgi_header, barcode5, barcode7, instrument, run):
    """
    Convert an MGI/BGI FASTQ header into an Illumina (Casava 1.8+) header.

    MGI headers have the form "<flowcell>L<lane>C<column>R<row><tile>/<read>",
    where <row> is three digits. The Illumina fields missing from an MGI
    header (instrument, run and index sequences) are supplied by the caller.
    The result has the form:

        @<instrument>:<run>:<flowcell>:<lane>:<tile>:<column>:<row> <read>:N:0:<barcode7>+<barcode5>

    Headers that already look like Illumina headers are returned unchanged
    (with a leading '@'). Anything after "/<read>" in an MGI header is
    discarded.

    Args:
        mgi_header (bytes): MGI header, with or without leading '@'.
        barcode5 (str): The i5 index sequence or barcode string.
        barcode7 (str): The i7 index sequence or barcode string.
        instrument (str): The sequencing instrument identifier.
        run (str): The run identifier or run number.

    Returns:
        bytes: Illumina-style header (including leading '@').

    Raises:
        ValueError: If the header matches neither the MGI nor the Illumina format.
    """
    mgi_header = mgi_header.lstrip(b"@").strip()
    if re.search(rb"^[\w-]+:\d+:[\w-]+:\d+:\d+:\d+:\d+\s+[12]:[YN]:\d+:", mgi_header):
        return b"@" + mgi_header
    strings = re.search(rb"^(\w+)L(\d+)C(\d+)R(\d{3})(\d+)\/([12])", mgi_header)
    if strings is None:
        raise ValueError(f"Header does not match expected MGI format: {mgi_header!r}")
    illumina_header = b"@%b:%b:%b:%d:%d:%d:%d %d:N:0:%b+%b" % (
        instrument.encode('ascii'),
        run.encode('ascii'),
        strings.group(1),
        int(strings.group(2)),
        int(strings.group(5)),
        int(strings.group(3)),
        int(strings.group(4)),
        int(strings.group(6)),
        barcode7.encode('ascii'),
        barcode5.encode('ascii')
    )
    return illumina_header
