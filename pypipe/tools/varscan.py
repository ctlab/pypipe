import sys

from pypipe import formats
from pypipe.utils import create_program, install_program


install_program("VarScan.sh", "VarScan.jar")


def pileup2snp(in_, out, min_coverage=None, min_reads2=None, log=None,
               min_avg_qual=None, min_var_freq=None, p_value=None):
    program = create_program("java -jar VarScan.jar pileup2snp",
            log, out, type_="jar")
    program.add_arg(in_, formats.Pileup)
    program.add_arg(min_coverage, int, "--min-coverage")
    program.add_arg(min_reads2, int, "--min-reads2")
    program.add_arg(min_avg_qual, int, "--min-avg-qual")
    program.add_arg(min_var_freq, float, "--min-var-freq")
    program.add_arg(p_value, float, "--p-value")
    return formats.Snp(out, program)


def pileup2indel(in_, out, min_coverage=None, min_reads2=None, log=None,
                 min_avg_qual=None, min_var_freq=None, p_value=None):
    program = create_program("java -jar VarScan.jar pileup2indel", log,
            out, type_="jar")
    program.add_arg(in_, formats.Pileup)
    program.add_arg(min_coverage, int, "--min-coverage")
    program.add_arg(min_reads2, int, "--min-reads2")
    program.add_arg(min_avg_qual, int, "--min-avg-qual")
    program.add_arg(min_var_freq, float, "--min-var-freq")
    program.add_arg(p_value, float, "--p-value")
    return formats.Indel(out, program)


def pileup2cns(in_, out, min_coverage=None, min_reads2=None, log=None,
               min_avg_qual=None, min_var_freq=None, p_value=None):
    program = create_program("java -jar VarScan.jar pileup2cns",
            log, out, type_="jar")
    program.add_arg(in_, formats.Pileup)
    program.add_arg(min_coverage, int, "--min-coverage")
    program.add_arg(min_reads2, int, "--min-reads2")
    program.add_arg(min_avg_qual, int, "--min-avg-qual")
    program.add_arg(min_var_freq, float, "--min-var-freq")
    program.add_arg(p_value, float, "--p-value")
    return formats.Cns(out, program)


def mpileup2snp(in_, out, min_coverage=None, min_reads2=None, log=None,
                min_avg_qual=None, min_var_freq=None, mon_freq_for_hom=None,
                p_value=None, strand_filter=None, output_vcf=None,
                variants=None):
    program = create_program("java -jar VarScan.jar mpileup2snp",
            log, out, type_="jar")
    program.add_arg(in_, formats.Pileup)
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
        return formats.Vcf(out, program)
    else:
        return formats.Snp(out, program)


def mpileup2indel(in_, out, min_coverage=None, min_reads2=None, log=None,
                  min_avg_qual=None, min_var_freq=None, mon_freq_for_hom=None,
                  p_value=None, strand_filter=None, output_vcf=None,
                  variants=None):
    program = create_program("java -jar VarScan.jar mpileup2indel",
            log, out, type_="jar")
    program.add_arg(in_, formats.Pileup)
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
        return formats.Vcf(out, program)
    else:
        return formats.Indel(out, program)


def mpileup2cns(in_, out, min_coverage=None, min_reads2=None, log=None,
                min_avg_qual=None, min_var_freq=None, mon_freq_for_hom=None,
                p_value=None, strand_filter=None, output_vcf=None,
                variants=None):
    program = create_program("java -jar VarScan.jar mpileup2cns",
            log, out, type_="jar")
    program.add_arg(_in, formats.Pileup)
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
        return formats.Vcf(out, program)
    else:
        return formats.Cns(out, program)


def somatic(normal, tumor, output, output_snp=None, output_indel=None,
            min_coverage=None, min_coverage_normal=None, log=None,
            min_coverage_tumor=None, min_var_freq=None, min_freq_for_hom=None,
            normal_purity=None, tumor_purity=None, p_value=None,
            somatic_p_value=None, strand_filter=None, validation=None):
    if output_snp and output_indel:
        sys.exit("Use only 'output_snp' or 'output_indel'")
    program = create_program("java -jar VarScan.jar somatic", log, None, "jar")
    program.add_arg(normal, formats.Pileup)
    program.add_arg(tumor, formats.Pileup)
    program.add_arg(output, str)
    program.add_arg(output_snp, bool, "--output-snp")
    program.add_arg(output_indel, bool, "--output-indel")
    program.add_arg(min_coverage, int, "--min-coverage")
    program.add_arg(min_coverage_normal, int, "--min-coverage-normal")
    program.add_arg(min_coverage_tumor, int, "--min-coverage-tumor")
    program.add_arg(min_var_freq, float, "--min-var-freq")
    program.add_arg(min_freq_for_hom, float, "--min-freq-for-hom")
    program.add_arg(normal_purity, float, "--normal-purity")
    program.add_arg(tumor_purity, float, "--tumor-purity")
    program.add_arg(p_value, float, "--p-value")
    program.add_arg(somatic_p_value, float, "--somatic-p-value")
    program.add_arg(strand_filter, bool, "--strand-filter")
    program.add_arg(validation, bool, "--validation")
    indel = formats.Indel(output + ".indel", program)
    snp = formats.Snp(output + ".snp", program)
    if output_snp:
        return indel
    elif output_indel:
        return snp
    return indel, snp


# This method doesnt work :(
def copynumber(normal, tumor, output, min_base_qual=None, min_map_qual=None,
               min_coverage=None, min_segment_size=None, max_segment_size=None,
               p_value=None, data_ratio=None, log=None):
    program = create_program("java -jar VarScan.jar copynumber", log,
            None, "jar")
    program.add_arg(normal, formats.Pileup)
    program.add_arg(tumor, formats.Pileup)
    program.add_arg(output, str)
    program.add_arg(min_base_qual, int, "--min-base-qual")
    program.add_arg(min_map_qual, int, "--min-map-qual")
    program.add_arg(min_coverage, int, "--min-coverage")
    program.add_arg(min_segment_size, int, "--min-segment-size")
    program.add_arg(max_segment_size, int, "--max-segment-size")
    program.add_arg(p_value, float, "--p-value")
    program.add_arg(data_ratio, float, "--data-ratio")
    return None 


# ... :(
