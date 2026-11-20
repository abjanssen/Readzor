"""Per-chunk processing that runs in the worker processes."""

import signal

from isal import igzip
import numpy as np

from readzor import constants
from readzor import stats as S
from readzor.constants import PHRED33_TO_64, PHRED64_TO_33
from readzor.arrays import average_quality_batch, qual_to_array, seq_to_array
from readzor.fastq_io import parse_blob, validate_fastq
from readzor.modules import adapter_trimming, adapter_trimming_overlap, average_quality_filter_wrapper, cut_set_ends, header_mgi_to_illumina, homopolymer_nucleotide_trimming, kmer_complexity_scan, sliding_window_quality, trim_ends_quality

def build_pipeline(parameters, read_direction=None):
    """
    Build the list of enabled trimming and filtering modules.

    Each entry is a (name, callable) pair. The name labels the module in the
    --stats output. The callable takes (sequence_arr, quality_arr,
    chunk_padding_bool, row_tilde_count, padding_mask_bool) for one chunk
    and returns a (left, right) pair of per-read boundary arrays, describing
    the [left, right) window of bases that module would keep (adapter
    trimming returns a third element with adapter hits for --stats). Callers merge
    all windows by taking the largest left and smallest right boundary per
    read, so the order of the modules does not affect the result.

    Modules included, when enabled in ``parameters`` (name in brackets):
      - k-mer complexity filter ("kmer_filter_flag") [kmer]
      - homopolymer end trimming ("poly_filter_flag") [poly]
      - adapter trimming ("adapter_filter_flag") [adapter]
      - quality-dependent end trimming ("endqual_filter_flag") [endqual]
      - input average-quality filter ("minimum_average_qual_pre" > 0) [avgqual_pre]
      - sliding-window quality trimming ("slider_filter_flag") [slider]
      - set-length end trimming ("cut_flag") [cut]
    Overlap trimming for paired reads is not part of this list; it is
    computed per pair and added separately in finish_reads().

    Args:
        parameters (dict): Run parameters (module flags and their settings).
        read_direction (str | None): "read_1" or "read_2" for one mate of a
            read pair, passed on to adapter trimming so that only that
            mate's directional adapters are searched for (see
            select_adapters_for_read()). None for unpaired reads.
            Defaults to None.

    Returns:
        list[tuple[str, callable]]: The enabled modules; empty if none are
            enabled.
    """
    pipeline = []
    #sequence_based
    if parameters.get("kmer_filter_flag"):
        pipeline.append(("kmer", lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: kmer_complexity_scan(seq, chunk_padding_bool, padding_mask_bool, kmer=parameters["kmer_size"], low_complex_cutoff=parameters["kmer_cutoff"], allow_n=parameters["allow_n_kmer"])))
    if parameters.get("poly_filter_flag"):
        pipeline.append(("poly", lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: homopolymer_nucleotide_trimming(seq, padding_mask_bool, row_tilde_count, chunk_padding_bool, poly_length_start=parameters["poly_length_start"], poly_length_end=parameters["poly_length_end"], poly_bases_start=parameters["poly_bases_start"], poly_bases_end=parameters["poly_bases_end"])))
    if parameters.get("adapter_filter_flag"): 
        pipeline.append(("adapter", lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: adapter_trimming(seq, chunk_padding_bool, row_tilde_count, adapter_sequences=parameters["adapter_sequences"], mismatches=parameters["adapter_mismatch"], adapter_seed=parameters["adapter_seed"], adapter_per_read = parameters["adapter_per_read"], read_direction=read_direction, adapter_group=None if parameters.get("adapter_excl") else parameters.get("adapter_group"))))
        
    #quality_based
    if parameters.get("endqual_filter_flag"):
        pipeline.append(("endqual", lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: trim_ends_quality(qual, chunk_padding_bool, row_tilde_count, padding_mask_bool, endqual_min_start=parameters["endqual_min_start"], endqual_min_end=parameters["endqual_min_end"])))
    if parameters.get("minimum_average_qual_pre") > 0:
        pipeline.append(("avgqual_pre", lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: average_quality_filter_wrapper(qual, chunk_padding_bool, row_tilde_count, min_avg_qual=parameters["minimum_average_qual_pre"])))
    if parameters.get("slider_filter_flag"):
        pipeline.append(("slider", lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: sliding_window_quality(qual, chunk_padding_bool, padding_mask_bool, row_tilde_count, slider_quality=parameters["slider_quality"], slider_window=parameters["slider_window"], slider_step=parameters["slider_step"])))

    #length_based
    if parameters.get("cut_flag"):
        pipeline.append(("cut", lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: cut_set_ends(seq, chunk_padding_bool, row_tilde_count, cut_start=parameters["cut_start"], cut_end=parameters["cut_end"])))

    return pipeline

def process_unpaired_chunk(chunk, phred_offset, minimum_average_qual_post, gzip_output, gzip_level, min_length_output, max_length_output, min_length_output_perc, max_length_output_perc, write_rejected, parameters, stats=None):
    """
    Validate, trim, filter and format one chunk of unpaired FASTQ reads.

    Steps:
      1. Each record is checked with validate_fastq(); invalid records are
         rejected.
      2. If MGI conversion is enabled, headers are converted to Illumina
         format and plus lines are reset to "+".
      3. All enabled modules from build_pipeline() are run. Their windows
         are merged by taking the largest left and smallest right boundary
         per read (the most stringent result per end).
      4. A read is kept if its merged window is non-empty, lies within the
         output length limits, and (if set) has a mean quality of at least
         minimum_average_qual_post.
      5. Kept reads are cut to their window and, if requested, their
         quality strings are re-encoded to parameters["phred_out"].

    Rejected records are collected untrimmed, with their original quality
    encoding. Records that failed validation are kept exactly as read;
    records rejected after trimming carry the converted header and plus
    line if MGI conversion is enabled.

    Args:
        chunk (Iterable[tuple[bytes, bytes, bytes, bytes]]): (header,
            sequence, plus, quality) records, as produced by parse_blob().
        phred_offset (int): Phred encoding offset (33 or 64).
        minimum_average_qual_post (float): Minimum acceptable mean quality
            after trimming. 0 disables this filter.
        gzip_output (bool): If True, compresses output records with isal gzip.
        gzip_level (int): isal gzip compression level (0-3).
        min_length_output (int | None): Minimum acceptable read length
            after trimming, as an absolute count.
        max_length_output (int | None): Maximum acceptable read length
            after trimming, as an absolute count.
        min_length_output_perc (int | None): Minimum acceptable read length
            after trimming, as a percentage of the read's input length.
            Ignored if min_length_output is given.
        max_length_output_perc (int | None): Maximum acceptable read length
            after trimming, as a percentage of the read's input length.
            Ignored if max_length_output is given.
        write_rejected (bool): If True, also collect rejected records
            (invalid on input or filtered out after trimming).
        parameters (dict): Run parameters, used to build the trimming
            pipeline and for N-filtering, input length limits,
            MGI-to-Illumina header conversion and Phred re-encoding.
        stats (dict | None): Read statistics for this chunk (--stats, see
            readzor.stats), filled in place. None to skip collection.

    Returns:
        tuple[bytes, int, int, bytes | list]: A tuple containing:
            - Kept records as FASTQ text (gzip-compressed if gzip_output).
            - Count of kept reads.
            - Count of rejected reads (invalid input records plus records
              filtered out after trimming).
            - Rejected records as FASTQ text (gzip-compressed if
              gzip_output) if write_rejected is True, otherwise an empty list.

    Raises:
        RuntimeError: If gzip compression fails.
    """
    valid_headers = []
    valid_sequences = []
    valid_pluses = []
    valid_qualities = []
    rejected_reads = []
    rejected = 0
    for header, sequence, plus, quality in chunk:
        if validate_fastq(header, sequence, plus, quality, n_filter = parameters["n_filter"], min_length_input = parameters["min_length_input"], max_length_input = parameters["max_length_input"]):
            valid_headers.append(header)
            valid_sequences.append(sequence)
            valid_pluses.append(plus)
            valid_qualities.append(quality)
        else:
            rejected += 1
            if write_rejected:
                rejected_reads.append(b"\n".join((header, sequence, plus, quality)) + b"\n")
    if stats is not None:
        S.add_invalid(stats, rejected)
    if not valid_headers:
        rejected_data = b"".join(rejected_reads) if write_rejected else []
        if gzip_output and rejected_data:
            rejected_data = igzip.compress(rejected_data, compresslevel=gzip_level)
        return b"", 0, rejected, rejected_data
    if parameters["mgi_convert_flag"]:
        valid_pluses = [b"+"] * len(valid_headers)
        valid_headers = [header_mgi_to_illumina(header, parameters["mgi_bc5"], parameters["mgi_bc7"], parameters["mgi_instrument"], parameters["mgi_run"]) for header in valid_headers]
    sequence_arr, chunk_padding_bool, row_tilde_count, padding_mask_bool, raw_lengths, max_len, n_reads = seq_to_array(sequence_list = valid_sequences)
    quality_arr = qual_to_array(quality_list = valid_qualities, phred_offset = phred_offset, chunk_padding_bool = chunk_padding_bool, padding_mask_bool = padding_mask_bool, max_len = max_len, n_reads = n_reads)
    n_reads, length = sequence_arr.shape
    left_list = [np.zeros(n_reads, dtype=np.int8)]
    right_list = [np.full(n_reads, length, dtype=np.int32) - row_tilde_count]
    module_names = []
    for name, step in build_pipeline(parameters):
        out = step(sequence_arr, quality_arr, chunk_padding_bool, row_tilde_count, padding_mask_bool)
        left_list.append(out[0])
        right_list.append(out[1])
        module_names.append(name)
        if stats is not None and len(out) > 2:
            S.collect_adapter_hits(stats, out[2]["adapters"], out[2]["hit_index"])
    lefts = np.maximum.reduce(left_list)
    rights = np.minimum.reduce(right_list)
    if (min_length_output is not None or max_length_output is not None
            or min_length_output_perc is not None or max_length_output_perc is not None):
        lengths_out = rights - lefts
        if max_length_output is not None:
            effective_max = max_length_output
        elif max_length_output_perc is not None:
            effective_max = (raw_lengths * max_length_output_perc / 100).astype(int)
        else:
            effective_max = np.inf
        if min_length_output is not None:
            effective_min = min_length_output
        elif min_length_output_perc is not None:
            effective_min = (raw_lengths * min_length_output_perc / 100).astype(int)
        else:
            effective_min = 0
        length_mask = (lengths_out <= effective_max) & (lengths_out >= effective_min) & (lengths_out > 0)
    else:
        lengths_out = rights - lefts
        length_mask = lengths_out > 0
    if minimum_average_qual_post > 0:
        average_quals = average_quality_batch(quality_arr, lefts, rights)
        qual_mask = average_quals >= minimum_average_qual_post
    else:
        qual_mask = np.ones(n_reads, dtype=bool)
    keep_mask = length_mask & qual_mask
    if stats is not None:
        S.collect_input(stats, quality_arr, raw_lengths)
        S.collect_modules(stats, module_names, left_list[1:], right_list[1:], raw_lengths, lefts, rights)
        S.collect_output(stats, quality_arr, lefts, rights, keep_mask)
    if parameters["phred_out"] == 33 and phred_offset == 64:
        qual_table = PHRED64_TO_33
    elif parameters["phred_out"] == 64 and phred_offset == 33:
        qual_table = PHRED33_TO_64
    else:
        qual_table = None
    results = []
    for header, sequence, plus_line, quality, left, right, keep in zip(
            valid_headers, valid_sequences, valid_pluses, valid_qualities, lefts, rights, keep_mask):
        if not keep:
            rejected += 1
            if write_rejected:
                rejected_reads.append(b"\n".join((header, sequence, plus_line, quality)) + b"\n")
            continue
        qual_out = quality[left:right]
        if qual_table is not None:
            qual_out = qual_out.translate(qual_table)
        results.append(b"\n".join((header, sequence[left:right], plus_line, qual_out)) + b"\n")
    len_results = len(results)
    results = b"".join(results)
    if write_rejected:
        rejected_reads = b"".join(rejected_reads)
    if gzip_output:
        try:
            results = igzip.compress(results, compresslevel=gzip_level)
            if write_rejected:
                rejected_reads = igzip.compress(rejected_reads, compresslevel=gzip_level)
        except Exception as e:
            raise RuntimeError(f"Compression failed inside worker: {str(e)}") from None
    return results, len_results, rejected, rejected_reads

def prepare_reads(records, phred_offset, write_rejected, parameters, stats=None):
    """
    Validate one mate's batch of FASTQ records and build the arrays every
    later step needs. Runs once per mate per chunk.

    Invalid records (see validate_fastq()) are rejected. If MGI conversion
    is enabled, the valid records' headers are converted to Illumina format
    and their plus lines reset to "+".

    Args:
        records (Iterable[tuple[bytes, bytes, bytes, bytes]]): (header,
            sequence, plus, quality) records, as produced by parse_blob().
        phred_offset (int): Phred encoding offset (33 or 64).
        write_rejected (bool): If True, also collect invalid records.
        parameters (dict): Run parameters (N-filtering, input length
            limits, MGI-to-Illumina header conversion).
        stats (dict | None): Read statistics of this mate (--stats); the
            invalid-record count is added to it. None to skip.

    Returns:
        tuple[dict | None, int, list[bytes]]:
            - The batch, or None if no record passed validation. The batch
              is a dict with:
                - "valid_reads" (numpy.ndarray): index in ``records`` of the
                  record in each array row.
                - "headers", "sequences", "pluses", "qualities" (list[bytes]):
                  the valid records' fields, one entry per array row.
                - "sequence_arr", "quality_arr", "chunk_padding_bool",
                  "row_tilde_count", "padding_mask_bool", "raw_lengths":
                  outputs of seq_to_array() / qual_to_array().
            - Count of records rejected during validation.
            - Rejected records, populated only if write_rejected is True.
    """
    valid_reads = []
    valid_headers = []
    valid_sequences = []
    valid_pluses = []
    valid_qualities = []
    rejected_reads = []
    rejected = 0
    for k, (header, sequence, plus, quality) in enumerate(records):
        if validate_fastq(header, sequence, plus, quality, n_filter = parameters["n_filter"], min_length_input = parameters["min_length_input"], max_length_input = parameters["max_length_input"]):
            valid_reads.append(k)
            valid_headers.append(header)
            valid_sequences.append(sequence)
            valid_pluses.append(plus)
            valid_qualities.append(quality)
        else:
            rejected += 1
            if write_rejected:
                rejected_reads.append(b"\n".join((header, sequence, plus, quality)) + b"\n")
    if stats is not None:
        S.add_invalid(stats, rejected)
    if not valid_headers:
        return None, rejected, rejected_reads
    if parameters["mgi_convert_flag"]:
        valid_pluses = [b"+"] * len(valid_headers)
        valid_headers = [header_mgi_to_illumina(header, parameters["mgi_bc5"], parameters["mgi_bc7"], parameters["mgi_instrument"], parameters["mgi_run"]) for header in valid_headers]
    sequence_arr, chunk_padding_bool, row_tilde_count, padding_mask_bool, raw_lengths, max_len, n_reads = seq_to_array(sequence_list = valid_sequences)
    quality_arr = qual_to_array(quality_list = valid_qualities, phred_offset = phred_offset, chunk_padding_bool = chunk_padding_bool, padding_mask_bool = padding_mask_bool, max_len = max_len, n_reads = n_reads)
    batch = {
        "valid_reads": np.array(valid_reads),
        "headers": valid_headers,
        "sequences": valid_sequences,
        "pluses": valid_pluses,
        "qualities": valid_qualities,
        "sequence_arr": sequence_arr,
        "quality_arr": quality_arr,
        "chunk_padding_bool": chunk_padding_bool,
        "row_tilde_count": row_tilde_count,
        "padding_mask_bool": padding_mask_bool,
        "raw_lengths": raw_lengths,
    }
    return batch, rejected, rejected_reads

def paired_overlap_cutoffs(batch_1, batch_2, parameters):
    """
    Run adapter_trimming_overlap() on the reads of a chunk that are valid in
    both mates, and map the results back onto each mate's own batch rows.

    Args:
        batch_1 (dict): Output of prepare_reads() for R1 of the chunk.
        batch_2 (dict): Output of prepare_reads() for R2 of the same chunk
            (record k of one is the mate of record k of the other). Neither
            batch may be None.
        parameters (dict): Run parameters; uses "overlap_min_length",
            "overlap_mismatches" and "overlap_portion_mismatch".

    Returns:
        tuple: ((left_1, right_1), (left_2, right_2)), one entry per batch
            row. Rows without a valid mate get left 0 and right int32 max,
            so the other modules decide their trimming alone.
    """
    no_cut = np.iinfo(np.int32).max
    n_1, n_2 = len(batch_1["headers"]), len(batch_2["headers"])
    left_1 = np.zeros(n_1, dtype=np.int32)
    right_1 = np.full(n_1, no_cut, dtype=np.int32)
    left_2 = np.zeros(n_2, dtype=np.int32)
    right_2 = np.full(n_2, no_cut, dtype=np.int32)
    _, rows_1, rows_2 = np.intersect1d(batch_1["valid_reads"], batch_2["valid_reads"], assume_unique=True, return_indices=True)
    if rows_1.size == 0:
        return (left_1, right_1), (left_2, right_2)
    seq_arr_1 = batch_1["sequence_arr"][rows_1]
    seq_arr_2 = batch_2["sequence_arr"][rows_2]
    pad_1, pad_2 = batch_1["chunk_padding_bool"], batch_2["chunk_padding_bool"]
    tilde_1 = batch_1["row_tilde_count"][rows_1] if pad_1 else 0
    tilde_2 = batch_2["row_tilde_count"][rows_2] if pad_2 else 0
    (l1, r1), (l2, r2) = adapter_trimming_overlap(seq_arr_1, seq_arr_2, pad_1, tilde_1, pad_2, tilde_2, min_overlap = parameters["overlap_min_length"], max_mismatch = parameters["overlap_mismatches"], max_mismatch_frac = parameters["overlap_portion_mismatch"])
    left_1[rows_1], right_1[rows_1] = l1, r1
    left_2[rows_2], right_2[rows_2] = l2, r2
    return (left_1, right_1), (left_2, right_2)

def finish_reads(batch, phred_offset, minimum_average_qual_post, min_length_output, max_length_output, min_length_output_perc, max_length_output_perc, write_rejected, overlap_cutoffs, parameters, read_direction=None, chunk_length = 0, stats=None):
    """
    Trim and length/quality-filter one prepared batch (one mate of a paired
    chunk) and format the survivors, placed at their record position in the
    chunk so they can be matched with the other mate's survivors.

    Runs all enabled modules from build_pipeline(), plus the overlap
    cutoffs if given, merges their windows (largest left, smallest right
    boundary), and applies the same keep rules as process_unpaired_chunk():
    non-empty window, within the output length limits, and (if set) mean
    quality of at least minimum_average_qual_post. Kept reads are cut to
    their window and optionally Phred re-encoded; rejected records are
    collected untrimmed.

    Args:
        batch (dict | None): Output of prepare_reads (None = nothing valid).
        phred_offset (int): Phred encoding offset (33 or 64).
        minimum_average_qual_post (float): Minimum acceptable mean quality
            after trimming. 0 disables this filter.
        min_length_output, max_length_output (int | None): Absolute length
            limits after trimming.
        min_length_output_perc, max_length_output_perc (int | None): Length
            limits as a percentage of input length; ignored if the absolute
            limit is given.
        write_rejected (bool): If True, also collect filtered-out records.
        overlap_cutoffs (tuple | None): (left, right) per batch row from
            paired_overlap_cutoffs, or None.
        parameters (dict): Run parameters (pipeline, Phred re-encoding).
        read_direction (str | None): "read_1" or "read_2", passed on to
            build_pipeline() for mate-specific adapter selection.
            Defaults to None.
        chunk_length (int): Number of records in the original chunk
            (valid and invalid). Defaults to 0.
        stats (dict | None): Read statistics of this mate (--stats), filled
            in place. None to skip collection. Defaults to None.

    Returns:
        tuple[list[bytes | None], int, list[bytes]]:
            - Survivors as a list of length ``chunk_length``: entry k is the
              formatted FASTQ record of chunk record k, or None if that
              record was invalid or filtered out.
            - Count of records filtered out here.
            - Rejected records (populated only if write_rejected is True).
    """
    rejected_reads = []
    rejected = 0
    survivors = [None] * chunk_length
    if batch is None:
        return survivors, rejected, rejected_reads
    original_index = batch["valid_reads"].tolist()
    valid_headers, valid_sequences = batch["headers"], batch["sequences"]
    valid_pluses, valid_qualities = batch["pluses"], batch["qualities"]
    sequence_arr, quality_arr = batch["sequence_arr"], batch["quality_arr"]
    chunk_padding_bool, row_tilde_count = batch["chunk_padding_bool"], batch["row_tilde_count"]
    padding_mask_bool, raw_lengths = batch["padding_mask_bool"], batch["raw_lengths"]
    n_reads, length = sequence_arr.shape
    left_list = [np.zeros(n_reads, dtype=np.int8)]
    right_list = [np.full(n_reads, length, dtype=np.int32) - row_tilde_count]
    module_names = []
    for name, step in build_pipeline(parameters, read_direction):
        out = step(sequence_arr, quality_arr, chunk_padding_bool, row_tilde_count, padding_mask_bool)
        left_list.append(out[0])
        right_list.append(out[1])
        module_names.append(name)
        if stats is not None and len(out) > 2:
            S.collect_adapter_hits(stats, out[2]["adapters"], out[2]["hit_index"])
    if overlap_cutoffs is not None:
        overlap_left, overlap_right = overlap_cutoffs
        left_list.append(overlap_left)
        right_list.append(overlap_right)
        module_names.append("overlap")
    lefts = np.maximum.reduce(left_list)
    rights = np.minimum.reduce(right_list)
    if (min_length_output is not None or max_length_output is not None
            or min_length_output_perc is not None or max_length_output_perc is not None):
        lengths_out = rights - lefts
    
        if max_length_output is not None:
            effective_max = max_length_output
        elif max_length_output_perc is not None:
            effective_max = (raw_lengths * max_length_output_perc / 100).astype(int)
        else:
            effective_max = np.inf
    
        if min_length_output is not None:
            effective_min = min_length_output
        elif min_length_output_perc is not None:
            effective_min = (raw_lengths * min_length_output_perc / 100).astype(int)
        else:
            effective_min = 0
        length_mask = (lengths_out <= effective_max) & (lengths_out >= effective_min) & (lengths_out > 0)
    else:
        lengths_out = rights - lefts
        length_mask = lengths_out > 0
    if minimum_average_qual_post > 0:
        avg_quals = average_quality_batch(quality_arr, lefts, rights)
        qual_mask = avg_quals >= minimum_average_qual_post
    else:
        qual_mask = np.ones(n_reads, dtype=bool)
    keep_mask = length_mask & qual_mask
    if stats is not None:
        S.collect_input(stats, quality_arr, raw_lengths)
        S.collect_modules(stats, module_names, left_list[1:], right_list[1:], raw_lengths, lefts, rights)
        S.collect_output(stats, quality_arr, lefts, rights, keep_mask)
    if parameters["phred_out"] == 33 and phred_offset == 64:
        qual_table = PHRED64_TO_33
    elif parameters["phred_out"] == 64 and phred_offset == 33:
        qual_table = PHRED33_TO_64
    else:
        qual_table = None
    for i, keep in enumerate(keep_mask):
        if not keep:
            rejected += 1
            if write_rejected:
                rejected_reads.append(b"\n".join((valid_headers[i], valid_sequences[i], valid_pluses[i], valid_qualities[i])) + b"\n")
            continue
        left, right = int(lefts[i]), int(rights[i])
        seq_out = valid_sequences[i][left:right]
        qual_out = valid_qualities[i][left:right]
        if qual_table is not None:
            qual_out = qual_out.translate(qual_table)
        survivors[original_index[i]] = b"\n".join((valid_headers[i], seq_out, valid_pluses[i], qual_out)) + b"\n"
    return survivors, rejected, rejected_reads

def process_paired_chunk(chunks, phred_offset_1, phred_offset_2, gzip_output, gzip_level, discard_singles, parameters, mate_stats=None):
    """
    Trim and filter one chunk of paired reads, then sort the survivors into
    intact pairs and orphaned singletons.

    Both mates are validated (prepare_reads()), optionally overlap-trimmed
    as a pair (paired_overlap_cutoffs(), when "overlap_filter_flag" is set),
    and then trimmed and filtered independently (finish_reads()). Survivors
    are matched by position in the chunk (record k of R1 with record k of
    R2): records whose mate also survived are output as pairs (in chunk
    order), the others as singletons.

    If parameters["interleaved_out"] is True (always the case with
    --stdout), pairs are written as alternating R1/R2 records into the first
    paired output, and both mates' rejected records are merged into the
    first rejected output; the second outputs are then empty.

    Args:
        chunks (tuple[Iterable, Iterable]): (chunk1, chunk2): R1 and R2
            records as (header, sequence, plus, quality) tuples (see
            parse_blob()), covering the same reads in the same order.
        phred_offset_1 (int): Phred encoding offset (33 or 64) for R1.
        phred_offset_2 (int): Phred encoding offset (33 or 64) for R2.
        gzip_output (bool): If True, compresses output records with isal gzip.
        gzip_level (int): isal gzip compression level (0-3).
        discard_singles (bool): If True, singleton records are dropped
            (their counts are still reported).
        parameters (dict): Run parameters, including "interleaved_out",
            "write_rejected", "overlap_filter_flag" and all module settings.
        mate_stats (dict | None): {"R1": stats, "R2": stats} read
            statistics (--stats), filled in place. None to skip collection.

    Returns:
        tuple: An 11-element tuple containing:
            - paired_out_1 (bytes): R1 records of surviving pairs (all pair
              records, alternating R1/R2, when interleaving).
            - paired_out_2 (bytes): R2 records of surviving pairs (empty when
              interleaving).
            - singles_out_1 (bytes): Surviving R1 records whose mate did not
              survive (empty if discard_singles).
            - singles_out_2 (bytes): Surviving R2 records whose mate did not
              survive (empty if discard_singles).
            - num_paired (int): Count of surviving read pairs.
            - num_R1_singles (int): Count of surviving R1 singletons.
            - num_R2_singles (int): Count of surviving R2 singletons.
            - rejected_1 (int): Count of rejected R1 reads (invalid plus
              filtered out).
            - rejected_2 (int): Count of rejected R2 reads.
            - rejected_R1 (bytes): Rejected R1 records, untrimmed (both
              mates' rejected records when interleaving); empty unless
              "write_rejected" is set.
            - rejected_R2 (bytes): Rejected R2 records (empty when
              interleaving).
        Non-empty byte outputs are gzip-compressed when gzip_output is True.

    Raises:
        ValueError: If the R1 and R2 chunks have different numbers of records.
        RuntimeError: If gzip compression fails.
    """
    chunk1, chunk2 = chunks
    chunk1, chunk2 = list(chunk1), list(chunk2)
    n_chunk = len(chunk1)
    if n_chunk != len(chunk2):
        raise ValueError("R1/R2 chunks have different numbers of reads.")
    stats_1 = mate_stats["R1"] if mate_stats is not None else None
    stats_2 = mate_stats["R2"] if mate_stats is not None else None
    batch_1, invalid_1, rejected_R1 = prepare_reads(chunk1, phred_offset_1, write_rejected = parameters["write_rejected"], parameters = parameters, stats = stats_1)
    batch_2, invalid_2, rejected_R2 = prepare_reads(chunk2, phred_offset_2, write_rejected = parameters["write_rejected"], parameters = parameters, stats = stats_2)
    overlap_1 = overlap_2 = None
    if parameters["overlap_filter_flag"] and batch_1 is not None and batch_2 is not None:
        overlap_1, overlap_2 = paired_overlap_cutoffs(batch_1, batch_2, parameters)
    survivors_1, filtered_1, filtered_R1 = finish_reads(batch_1, phred_offset_1, minimum_average_qual_post = parameters["minimum_average_qual_post"], min_length_output = parameters["min_length_output"], max_length_output = parameters["max_length_output"], min_length_output_perc = parameters["min_length_output_perc"], max_length_output_perc = parameters["max_length_output_perc"], write_rejected = parameters["write_rejected"], overlap_cutoffs = overlap_1, read_direction="read_1", parameters = parameters, chunk_length = n_chunk, stats = stats_1)
    survivors_2, filtered_2, filtered_R2 = finish_reads(batch_2, phred_offset_2, minimum_average_qual_post = parameters["minimum_average_qual_post"], min_length_output = parameters["min_length_output"], max_length_output = parameters["max_length_output"], min_length_output_perc = parameters["min_length_output_perc"], max_length_output_perc = parameters["max_length_output_perc"], write_rejected = parameters["write_rejected"], overlap_cutoffs = overlap_2, read_direction="read_2", parameters = parameters, chunk_length = n_chunk, stats = stats_2)
    rejected_1 = invalid_1 + filtered_1
    rejected_2 = invalid_2 + filtered_2
    rejected_R1 = rejected_R1 + filtered_R1
    rejected_R2 = rejected_R2 + filtered_R2

    paired_out_1 = []
    paired_out_2 = []
    singles_out_1 = []
    singles_out_2 = []
    for r1, r2 in zip(survivors_1, survivors_2):
        if r1 is not None:
            if r2 is not None:
                paired_out_1.append(r1)
                paired_out_2.append(r2)
            else:
                singles_out_1.append(r1)
        elif r2 is not None:
            singles_out_2.append(r2)
    num_R1_singles = len(singles_out_1)
    num_R2_singles = len(singles_out_2)
    num_paired = len(paired_out_1)
    
    if discard_singles:
        singles_out_1 = b""
        singles_out_2 = b""
        
    interleave = parameters["interleaved_out"]
    if interleave:
        paired_out_1 = b"".join(r1 + r2 for r1, r2 in zip(paired_out_1, paired_out_2))
        paired_out_2 = b""
        rejected_R1 = b"".join(rejected_R1) + b"".join(rejected_R2)
        rejected_R2 = b""
    else:
        paired_out_1 = b"".join(paired_out_1)
        paired_out_2 = b"".join(paired_out_2)
        rejected_R1 = b"".join(rejected_R1)
        rejected_R2 = b"".join(rejected_R2)
    singles_out_1 = b"".join(singles_out_1)
    singles_out_2 = b"".join(singles_out_2)

    if gzip_output:
        try:
            if paired_out_1:
                paired_out_1 = igzip.compress(paired_out_1, compresslevel=gzip_level)
            if paired_out_2:
                paired_out_2 = igzip.compress(paired_out_2, compresslevel=gzip_level)
            if rejected_R1:
                rejected_R1 = igzip.compress(rejected_R1, compresslevel=gzip_level)
            if rejected_R2:
                rejected_R2 = igzip.compress(rejected_R2, compresslevel=gzip_level)
            if singles_out_1:
                singles_out_1 = igzip.compress(singles_out_1, compresslevel=gzip_level)
            if singles_out_2:
                singles_out_2 = igzip.compress(singles_out_2, compresslevel=gzip_level)
        except Exception as e:
            raise RuntimeError(f"Compression failed inside worker: {str(e)}") from None
    return paired_out_1, paired_out_2, singles_out_1, singles_out_2, num_paired, num_R1_singles, num_R2_singles, rejected_1, rejected_2, rejected_R1, rejected_R2

def process_unpaired_task_flat(task, parameters, mate_stats=None):
    """
    Worker wrapper for an unpaired task: parses the task's blob
    (parse_blob()), runs process_unpaired_chunk() on it and tags the result
    with the task type and filepath.

    Args:
        task (dict): An "unpaired" task from generate_unpaired_tasks().
        parameters (dict): Run parameters.
        mate_stats (dict | None): {"all": stats} read statistics (--stats),
            or None.

    Returns:
        tuple: (type, filepath, chunk_results, kept, rejected, rejected_reads),
            where chunk_results and rejected_reads are the kept and rejected
            record data, and kept and rejected are counts (see
            process_unpaired_chunk()).
    """
    chunk_results, kept, rejected, rejected_reads = process_unpaired_chunk(
        chunk=parse_blob(task["chunk"]),
        phred_offset=task["phred_offset"],
        minimum_average_qual_post=parameters["minimum_average_qual_post"],
        gzip_output = parameters["gzip_output"],
        gzip_level = parameters["gzip_level"],
        min_length_output = parameters["min_length_output"],
        max_length_output = parameters["max_length_output"],
        min_length_output_perc = parameters["min_length_output_perc"],
        max_length_output_perc =parameters["max_length_output_perc"],
        write_rejected= parameters["write_rejected"],
        parameters=parameters,
        stats=mate_stats["all"] if mate_stats is not None else None
    )
    return task["type"], task["filepath"], chunk_results, kept, rejected, rejected_reads

def process_paired_task_flat(task, parameters, mate_stats=None):
    """
    Worker wrapper for a paired task (two-file or interleaved): parses the
    task's blob(s) (parse_blob()), runs process_paired_chunk() on the R1/R2
    records and tags the result with the task type and filepaths. For an
    interleaved task, the records of its single blob are split alternately
    into R1 (records 1, 3, 5, ...) and R2 (records 2, 4, 6, ...).

    Args:
        task (dict): A "paired" task from generate_paired_tasks() or
            generate_unpaired_tasks(filetype="interleaved").
        parameters (dict): Run parameters.
        mate_stats (dict | None): {"R1": stats, "R2": stats} read
            statistics (--stats), or None.

    Returns:
        tuple: (type, file1, file2, paired_out_1, paired_out_2,
            R1_singles_out, R2_singles_out, num_paired, num_R1_singles,
            num_R2_singles, rejected_1, rejected_2, rejected_R1, rejected_R2).
            rejected_1/rejected_2 are rejected-read counts and
            rejected_R1/rejected_R2 the rejected record data, kept separate
            per mate since each mate is trimmed and filtered independently
            before reconciliation (see process_paired_chunk()).
    """
    file1 = task["file1"]
    file2 = task["file2"]
    if task.get("interleaved"):
        records = list(parse_blob(task["chunk1"]))
        chunk1, chunk2 = records[0::2], records[1::2]
    else:
        chunk1, chunk2 = parse_blob(task["chunk1"]), parse_blob(task["chunk2"])
    paired_out_1, paired_out_2, R1_singles_out, R2_singles_out, num_paired, num_R1_singles, num_R2_singles, rejected_1, rejected_2, rejected_R1, rejected_R2 = process_paired_chunk(
        chunks = (chunk1, chunk2),
        phred_offset_1=task["phred_offset_1"],
        phred_offset_2=task["phred_offset_2"],
        gzip_output = task["gzip_output"],
        gzip_level = task["gzip_level"],
        discard_singles=task["discard_singles"],
        parameters=parameters,
        mate_stats=mate_stats
    )
    return task["type"], file1, file2, paired_out_1, paired_out_2, R1_singles_out, R2_singles_out, num_paired, num_R1_singles, num_R2_singles, rejected_1, rejected_2, rejected_R1, rejected_R2

def worker_initializer(parameters):
    """
    Pool initializer: runs once in each worker process at startup.

    Makes the worker ignore SIGINT, so Ctrl+C is handled by the main
    process only, and stores `parameters` in the WORKER_PARAMETERS global,
    so it is pickled and sent once per worker instead of with every task.

    Args:
        parameters (dict): Run parameters.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    constants.WORKER_PARAMETERS = parameters

def unified_worker(task):
    """
    Pool entry point: route a task to the paired or unpaired handler, using
    the run parameters stored by worker_initializer().

    Args:
        task (dict): A task from generate_unpaired_tasks() or
            generate_paired_tasks(). Its "type" field ("paired" or
            "unpaired") selects the handler.

    Returns:
        tuple: The result of process_paired_task_flat() or
            process_unpaired_task_flat(), with one more element appended:
            the chunk's read statistics (see readzor.stats.new_chunk_report())
            with --stats, otherwise None.

    Raises:
        ValueError: If the task type is unknown or missing.
    """
    task_type = task.get("type")
    parameters = constants.WORKER_PARAMETERS
    report = None
    if task_type == "paired":
        if parameters["stats"]:
            report = S.new_chunk_report(("R1", "R2"))
        result = process_paired_task_flat(task, parameters=parameters, mate_stats=report["mates"] if report else None)
    elif task_type == "unpaired":
        if parameters["stats"]:
            report = S.new_chunk_report(("all",))
        result = process_unpaired_task_flat(task, parameters=parameters, mate_stats=report["mates"] if report else None)
    else:
        raise ValueError(f"Unknown or missing task type: {task_type}")
    return result + (report,)
