#!/usr/bin/env python3

from argparse import ArgumentParser
from pythia.models import load_model
from pythia.utils import predict_across
import pandas as pd
from loguru import logger
import tensorflow as tf
import os

## Set GPU device. Change number for different jobs.
## Otherwise tf automatically uses the lowest numbered GPU

CHUNKSIZE = 1000000

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("model")
    parser.add_argument("tiled")
    parser.add_argument("outfile")

    args = parser.parse_args()
    
    ## Set GPU device. Change number for different jobs.
    ## Otherwise tf automatically uses the lowest numbered GPU
    os.environ["CUDA_VISIBLE_DEVICES"]="3"
    with tf.device('/device:GPU:3'):
        model = load_model(args.model)

        predicted = []
        # chunk to handle huge files
        tiled = pd.read_csv(
            args.tiled, delimiter="\t", names=["sequence", "region"], chunksize=CHUNKSIZE
        )
        logger.info(
            f"Predicting activity of sequences from {args.tiled} using {args.model}."
        )
    
        nprocessed = 0
        for chunk in tiled:
            predicted_chunk = predict_across(model, chunk.sequence)
            predicted_chunk = predicted_chunk.merge(chunk, how="left", on="sequence")
            nprocessed += len(chunk)
            logger.info(f"Processed {nprocessed} sequences.")
            if nprocessed == CHUNKSIZE:
                predicted_chunk.to_csv(args.outfile, sep="\t", index=False)
            else:
                predicted_chunk.to_csv(
                    args.outfile, sep="\t", index=False, mode="a", header=False
                )
                
    logger.info(
        f"Finished predicting activity. Predicted activity values are in {args.outfile}."
    )
