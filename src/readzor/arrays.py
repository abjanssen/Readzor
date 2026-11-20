"""Conversion of reads into numpy arrays, and array helpers."""

import numpy as np

def seq_to_array(sequence_list):
    """
    Convert a list of nucleotide sequences into a 2D numpy array of ASCII
    codes, and work out the chunk's padding layout for all downstream steps.

    If all reads have the same length, the sequences are joined directly.
    Otherwise each read is right-padded with null bytes (0x00) to the length
    of the longest read. The raw bytes are then reinterpreted as a
    (n_reads, max_len) matrix for vectorized processing.

    Args:
        sequence_list (list[bytes]): Sequence byte-strings (non-empty list).

    Returns:
        tuple: A tuple containing:
            - numpy.ndarray: Signed 8-bit integer array of shape
              (n_reads, max_len) of ASCII character codes.
            - bool: chunk_padding_bool, True if reads differ in length.
            - numpy.ndarray | int: row_tilde_count, (n_reads,) number of
              padding bytes per read, or 0 when not padded.
            - numpy.ndarray | None: padding_mask_bool, (n_reads, max_len)
              True at padding positions, or None when not padded.
            - numpy.ndarray: (n_reads,) real read lengths.
            - int: max_len, length of the longest read.
            - int: n_reads, number of reads.
    """
    sequence_arr_lengths = np.array([len(s) for s in sequence_list])
    n_reads = len(sequence_list)
    max_len = sequence_arr_lengths.max()
    min_len = sequence_arr_lengths.min()
    if max_len == min_len:
        joined = b''.join(sequence_list)
        array = np.frombuffer(joined, dtype=np.int8).reshape(n_reads, max_len)
        chunk_padding_bool = False
        row_tilde_count = 0
        padding_mask_bool = None
        return array, chunk_padding_bool, row_tilde_count, padding_mask_bool, sequence_arr_lengths, max_len, n_reads
    else:
        padded = [s.ljust(max_len, b'\x00') for s in sequence_list]
        joined = b''.join(padded)
        array = np.frombuffer(joined, dtype=np.int8).reshape(n_reads, max_len)
        padding_mask_bool = array == ord('\x00')
        row_tilde_count = np.sum(padding_mask_bool, axis=1)
        chunk_padding_bool = bool(np.any(row_tilde_count > 0))
        return array, chunk_padding_bool, row_tilde_count, padding_mask_bool, sequence_arr_lengths, max_len, n_reads

def qual_to_array(quality_list, phred_offset, chunk_padding_bool, padding_mask_bool, max_len, n_reads):
    """
    Convert a list of Phred quality strings into a 2D numeric numpy array,
    using the padding layout seq_to_array() determined for the same reads.

    Concatenates all quality strings (right-padded with null bytes when
    chunk_padding_bool is True), reinterprets the raw bytes as integers,
    subtracts the Phred offset (padding positions stay 0), and reshapes
    into a (n_reads, max_len) matrix.

    Args:
        quality_list (list[bytes]): Quality byte-strings, each the same length
            as its sequence (guaranteed by validate_fastq).
        phred_offset (int): The Phred encoding offset (33 or 64) to subtract.
        chunk_padding_bool (bool): From seq_to_array.
        padding_mask_bool (numpy.ndarray | None): From seq_to_array.
        max_len (int): From seq_to_array.
        n_reads (int): From seq_to_array.

    Returns:
        numpy.ndarray: Signed 8-bit integer array of shape
            (n_reads, max_len) with true quality scores.
    """
    if not chunk_padding_bool:
        joined = b''.join(quality_list)
        array = np.frombuffer(joined, dtype=np.int8).reshape(n_reads, max_len)
        return array - phred_offset
    padded = [q.ljust(max_len, b'\x00') for q in quality_list]
    joined = b''.join(padded)
    array = np.frombuffer(joined, dtype=np.int8).reshape(n_reads, max_len)
    return np.where(padding_mask_bool, array, array - phred_offset)

def average_quality_batch(quality_arr, lefts, rights):
    """
    Compute the mean quality score within a per-read [left, right) window,
    vectorized across all reads at once.

    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base quality scores.
        lefts (numpy.ndarray | int): Per-read left boundary (inclusive),
            shape (n_reads,), or a single value for all reads.
        rights (numpy.ndarray | int): Per-read right boundary (exclusive),
            shape (n_reads,), or a single value for all reads.

    Returns:
        numpy.ndarray: Per-read mean quality within [left, right), shape
            (n_reads,). Reads with an empty window (right <= left) get 0.0.
    """
    n_reads, length = quality_arr.shape
    cumsum = np.empty((n_reads, length + 1), dtype=np.int32)
    cumsum[:, 0] = 0
    np.cumsum(quality_arr, axis=1, dtype=np.int32, out=cumsum[:, 1:])
    row_indices = np.arange(n_reads)
    sums = cumsum[row_indices, rights] - cumsum[row_indices, lefts]
    counts = rights - lefts
    avg_qualities = np.divide(sums, counts, out=np.zeros_like(sums, dtype=np.float64), where=counts > 0)
    return avg_qualities
