#!/usr/bin/env python3

from pysam import FastaFile
import toolz as tz
from loguru import logger
from argparse import ArgumentParser

WINDOW_SIZE = 50
SEQUENCE_SIZE = 170
VALID_CHROMOSOMES = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chrX",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr20",
    "chrY",
    "chr19",
    "chr22",
    "chr21",
]


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("reference", help="Path to indexed reference file.")

    args = parser.parse_args()

    reference = args.reference
    ff = FastaFile(reference)

    for chromosome in ff.references:
        #if chromosome not in VALID_CHROMOSOMES:
        #    logger.info(f"Skipping {chromosome} as it is not a main chromosome.")
        #    continue
        logger.info(f"Processing {chromosome}.")
        chrom = chromosome
        chromlen = ff.get_reference_length(chrom)
        start = 0
        end = SEQUENCE_SIZE
        processed = 0
        skipped = 0
        while end < chromlen:
            processed += 1
            region = f"{chrom}:{start}-{end}"
            if processed % 100000 == 0:
                logger.info(f"On {region}")
                logger.info(f"Processed {processed} windows.")
            sequence = ff.fetch(reference=chrom, start=start, end=end)
            if "N" in sequence:
                skipped += 1
                start = start + WINDOW_SIZE
                end = end + WINDOW_SIZE
                continue
            print(f"{sequence.upper()}\t{region}")
            start = start + WINDOW_SIZE
            end = end + WINDOW_SIZE
        logger.info(f"Skipped {skipped} windows in {chromosome} due to containing Ns.")
