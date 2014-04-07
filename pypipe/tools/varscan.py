import sys

from pypipe import formats
from pypipe.utils import create_program, install_program


install_program("VarScan.sh", "VarScan.jar")


def pileup2snp(_in, _out, min_coverage=None, min_reads2=None, log=None,
               min_avg_qual=None, min_var_freq=None, p_value=None):
    program = create_program("java -jar VarScan.jar pileup2snp",
            log, _out, type_="jar")
    program.add_args(_in, formats.Pileup, 1)
    program.add_arg(min_coverage, int, "--min-coverage")
    program.add_arg(min_reads2, int, "--min-reads2")
    program.add_arg(min_avg_qual, int, "--min-avg-qual")
    program.add_arg(min_var_freq, float, "--min-var-freq")
    program.add_arg(p_value, float, "--p-value")
    return formats.Snp(_out, program)


def pileup2indel(_in, _out, min_coverage=None, min_reads2=None, log=None,
                 min_avg_qual=None, min_var_freq=None, p_value=None):
    program = create_program("java -jar VarScan.jar pileup2indel", log,
            _out, type_="jar")
    program.add_args(_in, formats.Pileup, 1)
    program.add_arg(min_coverage, int, "--min-coverage")
    program.add_arg(min_reads2, int, "--min-reads2")
    program.add_arg(min_avg_qual, int, "--min-avg-qual")
    program.add_arg(min_var_freq, float, "--min-var-freq")
    program.add_arg(p_value, float, "--p-value")
    return formats.Indel(_out, program)


def pileup2cns(_in, _out, min_coverage=None, min_reads2=None, log=None,
               min_avg_qual=None, min_var_freq=None, p_value=None):
    program = create_program("java -jar VarScan.jar pileup2cns",
            log, _out, type_="jar")
    program.add_args(_in, formats.Pileup, 1)
    program.add_arg(min_coverage, int, "--min-coverage")
    program.add_arg(min_reads2, int, "--min-reads2")
    program.add_arg(min_avg_qual, int, "--min-avg-qual")
    program.add_arg(min_var_freq, float, "--min-var-freq")
    program.add_arg(p_value, float, "--p-value")
    return formats.Cns(_out, program)


def mpileup2snp(_in, _out, min_coverage=None, min_reads2=None, log=None,
                min_avg_qual=None, min_var_freq=None, mon_freq_for_hom=None,
                p_value=None, strand_filter=None, output_vcf=None,
                variants=None):
    program = create_program("java -jar VarScan.jar mpileup2snp",
            log, _out, type_="jar")
    program.add_args(_in, formats.Pileup, 1)
    program.add_arg(min_coverage, int, "--min-coverage")
    program.add_arg(min_reads2, int, "--min-reads2")
    program.add_arg(min_avg_qual, int, "--min-avg-qual")
    program.add_arg(min_var_freq, float, "--min-val-freq")
    program.add_arg(min_freq_for_hom, float, "--min-freq-for-hom")
    program.add_arg(p_value, float, "--p-value")
    program.add_arg(strand_filter, int, "--strand-filter")
    program.add_arg(output_vcf, bool, "--output-vcf")
    program.add_arg(variants, int, "--variants")
    if output_vcf:
        return formats.Vcf(_out, program)
    else:
        return formats.Snp(_out, program)


def mpileup2indel(_in, _out, min_coverage=None, min_reads2=None, log=None,
                  min_avg_qual=None, min_var_freq=None, mon_freq_for_hom=None,
                  p_value=None, strand_filter=None, output_vcf=None,
                  variants=None):
    program = create_program("java -jar VarScan.jar mpileup2indel",
            log, _out, type_="jar")
    program.add_args(_in, formats.Pileup, 1)
    program.add_arg(min_coverage, int, "--min-coverage")
    program.add_arg(min_reads2, int, "--min-reads2")
    program.add_arg(min_avg_qual, int, "--min-avg-qual")
    program.add_arg(min_var_freq, float, "--min-val-freq")
    program.add_arg(min_freq_for_hom, float, "--min-freq-for-hom")
    program.add_arg(p_value, float, "--p-value")
    program.add_arg(strand_filter, int, "--strand-filter")
    program.add_arg(output_vcf, bool, "--output-vcf")
    program.add_arg(variants, int, "--variants")
    if output_vcf:
        return formats.Vcf(_out, program)
    else:
        return formats.Indel(_out, program)


def mpileup2cns(_in, _out, min_coverage=None, min_reads2=None, log=None,
                min_avg_qual=None, min_var_freq=None, mon_freq_for_hom=None,
                p_value=None, strand_filter=None, output_vcf=None,
                variants=None):
    program = create_program("java -jar VarScan.jar mpileup2cns",
            log, _out, type_="jar")
    program.add_args(_in, formats.Pileup, 1)
    program.add_arg(min_coverage, int, "--min-coverage")
    program.add_arg(min_reads2, int, "--min-reads2")
    program.add_arg(min_avg_qual, int, "--min-avg-qual")
    program.add_arg(min_var_freq, float, "--min-val-freq")
    program.add_arg(min_freq_for_hom, float, "--min-freq-for-hom")
    program.add_arg(p_value, float, "--p-value")
    program.add_arg(strand_filter, int, "--strand-filter")
    program.add_arg(output_vcf, bool, "--output-vcf")
    program.add_arg(variants, int, "--variants")
    if output_vcf:
        return formats.Vcf(_out, program)
    else:
        return formats.Cns(_out, program)


def somatic(normal, tumor, output, output_snp=None, output_indel=None,
            min_coverage=None, min_coverage_normal=None,
            min_coverage_tumor=None, min_var_freq=None, min_freq_for_hom=None,
            normal_purity=None, tumor_purity=None, p_value=None,
            somatic_p_value=None, strand_filter=None, validation=None):
    # TODO need for speed
    return None

