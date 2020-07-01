#!/usr/bin/env python3

import argparse
from hfad_fitter import run as hfad_fit
from mpi4py import MPI, futures


def cmdargs():
    parser = argparse.ArgumentParser(
        description="Fit high frequency atmospheric distortion (hfad)")
    parser.add_argument("in_file_list_path", type=str, default=None,
                        help="Path to list of input file paths")
    parser.add_argument("out_file_list_path", type=str, default=None,
                        help="Path to list of output file paths")
    return parser.parse_args()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def run(in_file_path, out_file_path):
    try:
        res = hfad_fit(file_path, out_file_path)
    except:
        res = False
    return res


if __name__=="__main__":
    # grab command line arguments
    clargs = cmdargs()

    # grab file lists
    with open(clargs.in_file_list_path, "r") as f:
        infiles = [line.strip() for line in f.readlines()]
    with open(clargs.out_file_list_path, "r") as f:
        outfiles = [line.strip() for line in f.readlines()]

    # pair up input and output files
    in_out = list(zip(infiles, outfiles))

    # run the fitter asynchronously across all available nodes
    with futures.MPICommExecutor(comm=comm, root=0) as executor:
        if executor is not None:
            ret = list(executor.map(run, in_out))

    # check that all files processed correctly
    assert sum(ret)==len(files)
