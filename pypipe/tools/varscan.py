import sys

from pypipe import formats
from pypipe.utils import create_program, install_program


install_program("VarScan.sh", "VarScan.jar")


def pileup2snp(_in, _out, min_coverage=None, min_reads2=None, log=None,
               min_avg_qual=None, min_var_freq=None, p_value=None):
    program = create_program("java -jar VarScan.jar", log, _out)
    program.add_args(_in, formats.Pileup, 1)
    program.add_arg(min_coverage, int, "--min-coverage")
    program.add_arg(min_reads2, int, "--min-reads2")
    program.add_arg(min_avg_qual, int, "--min-avg-qual")
    program.add_arg(min_var_freq, float, "--min-val-freq")
    program.add_arg(p_value, float, "--p-value")
    return formats.Snp(_out, program)

