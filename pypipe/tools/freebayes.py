import sys

from pypipe import formats
from pypipe.utils import create_program, install_program


install_program("freebayes.sh", "freebayes")


def freebayes(_in, v, f, vcf=None, b=None, bam=None, L=None, bam_list=None,
              fasta_reference=None, t=None, targets=None, r=None,
              region=None, s=None, samples=None, populations=None, A=None,
              cnv_map=None, trace=None, failed_alleles=None,
              variant_input=None, l=None, only_use_input_alleles=None,
              haplotype_basis_alleles=None, report_all_haplotype_alleles=None,
              report_monorphic=None, P=None, pvar=None, T=None, theta=None,
              p=None, ploidy=None, J=None, pooled_discrete=None, K=None,
              pooled_continuous=None, Z=None, use_reference_allele=None,
              reference_quality=None, I=None, no_snps=None, i=None,
              no_indels=None, X=None, no_mnps=None, u=None, no_complex=None,
              n=None, use_best_n_alleles=None, E=None, max_complex_gap=None,
              haplotype_length=None, min_repeat_length=None,
              min_repeat_entropy=None, no_partial_observations=None, O=None,
              dont_left_align_indels=None, _4=None, use_duplicate_reads=None,
              m=None, min_mapping_quality=None, q=None, min_base_quality=None,
              R=None, min_supporting_allele_qsum=None, Y=None,
              min_supporting_mapping_qsum=None, Q=None,
              mismatch_base_quality_threshold=None, U=None,
              read_mismatch_limit=None, z=None, log=None,
              read_max_mismatch_fraction=None, read_snp_limit=None, e=None,
              read_indel_limit=None, _0=None, standard_filters=None,
              F=None, min_alternate_fraction=None, C=None,
              min_alternate_count=None, _3=None, min_alternate_qsum=None,
              G=None, min_alternate_total=None, min_coverage=None, k=None,
              no_population_priors=None, w=None, hwe_priors_off=None, V=None,
              binomial_obs_priors_off=None, a=None,
              allele_balance_priors_off=None, observation_bias=None,
              base_quality_cap=None, experimental_gls=None,
              prob_contamination=None, contamination_estimates=None,
              report_genotype_likelihood_max=None, B=None,
              genotyping_max_iterations=None, genotyping_max_banddepth=None,
              W=None, posterior_integration_limits=None, N=None,
              exclude_unobserved_genotypes=None, S=None,
              genotype_variant_threshold=None, j=None,
              use_mapping_quality=None, H=None, harmonic_indel_quality=None,
              D=None, read_dependence_factor=None, genotype_qualities=None):
    if v and vcf:
        sys.exit("Use only 'v' or 'vcf'")
    if b and bam:
        sys.exit("Use only 'b' or 'bam'")
    if L and bam_list:
        sys.exit("Use only 'L' or 'bam_list'")
    if f and fasta_reference:
        sys.exit("Use only 'f' or 'fasta_reference'")
    if t and targets:
        sys.exit("Use only 't' or 'targets'")
    if r and region:
        sys.exit("Use only 'r' or 'region'")
    if s and samples:
        sys.exit("Use only 's' or 'samples'")
    if A and cnv_map:
        sys.exit("Use only 'A' or 'cnv_map'")
    if l and only_use_input_alleles:
        sys.exit("Use only 'l' or 'only_use_input_alleles'")
    if P and pvar:
        sys.exit("Use only 'P' or 'pvar'")
    if T and theta:
        sys.exit("Use only 'T' or 'theta'")
    if p and ploidy:
        sys.exit("Use only 'p' or 'ploidy'")
    if J and pooled_discrete:
        sys.exit("Use only 'J' or 'pooled_discrete'")
    if K and pooled_continuous:
        sys.exit("Use only 'K' or 'pooled_continuous'")
    if Z and use_reference_allele:
        sys.exit("Use only 'Z' or 'use_reference_allele'")
    if I and no_snps:
        sys.exit("Use only 'I' or 'no_snps'")
    if i and no_indels:
        sys.exit("Use only 'i' or 'no_indels'")
    if X and no_mnps:
        sys.exit("Use only 'X' or 'no_mnps'")
    if u and no_complex:
        sys.exit("Use only 'u' or 'no_complex'")
    if n and use_best_n_alleles:
        sys.exit("Use only 'n' or 'use_best_n_alleles'")
    if E and max_complex_gap or E and haplotype_length or \
            max_complex_gap and haplotype_length:
        sys.exit("Use only 'E' or 'max_complex_gap' or 'haplotype_length'")
    if O and dont_left_align_indels:
        sys.exit("Use only 'O' or 'dont_left_align_indels'")
    if _4 and use_duplicate_reads:
        sys.exit("Use only '_4' or 'use_only_duplicate_reads'")
    if m and min_mapping_quality:
        sys.exit("Use only 'm' or 'min_mapping_quality'")
    if q and min_base_quality:
        sys.exit("Use only 'q' or 'min_base_quality'")
    if R and min_supporting_allele_qsum:
        sys.exit("Use only 'R' or 'min_supporting_allele_qsum'")
    if Y and min_supporting_mapping_qsum:
        sys.exit("Use only 'Y' or 'min_supporting_mapping_qsum'")
    if Q and mismatch_base_quality_threshold:
        sys.exit("Use only 'Q' or 'mismatch_base_quality_threshold'")
    if U and read_mismatch_limit:
        sys.exit("Use only 'U' or 'read_mismatch_limit'")
    if z and read_max_mismatch_fraction:
        sys.exit("Use only 'z' or 'read_max_mismatch_fraction'")
    if e and read_indel_limit:
        sys.exit("Use only 'e' or 'read_indel_limit'")
    if _0 and standard_filters:
        sys.exit("Use only '_0' or 'standard_filters'")
    if F and min_alternate_fraction:
        sys.exit("Use only 'F' or 'min_alternate_fraction'")
    if C and min_alternate_count:
        sys.exit("Use only 'C' or 'min_alternate_count'")
    if _3 and min_alternate_qsum:
        sys.exit("Use only '_3' or 'min_alternate_qsum'")
    if G and min_alternate_total:
        sys.exit("Use only 'G' or 'min_alternate_total'")
    if k and no_population_priors:
        sys.exit("Use only 'k' or 'no_population_priors'")
    if w and hwe_priors_off:
        sys.exit("Use only 'w' or 'hwe_priors_off'")
    if V and binomial_obs_priors_off:
        sys.exit("Use only 'V' or 'binomial_obs_priors_off'")
    if a and allele_balance_priors_off:
        sys.exit("Use only 'a' or 'allele_balance_priors_off'")
    if B and genotyping_max_banddepth:
        sys.exit("Use only 'B' or 'genotyping_max_banddepth'")
    if W and posterior_integration_limits:
        sys.exit("Use only 'W' or 'posterior_integration_limits'")
    if N and exclude_unobserved_genotypes:
        sys.exit("Use only 'N' or 'exclude_unobserved_genotypes'")
    if S and genotype_variant_threshold:
        sys.exit("Use only 'S' or 'genotype_variant_threshold'")
    if j and use_mapping_quality:
        sys.exit("Use only 'j' or 'use_mapping_quality'")
    if H and harmonic_indel_quality:
        sys.exit("Use only 'H' or 'harmonic_indel_quality'")
    if D and read_dependence_factor:
        sys.exit("Use only 'D' or 'read_dependence_factor'")
    program = create_program("freebayes", log)
    program.add_arg(b, formats.Bam, "-b")
    program.add_arg(bam, formats.Bam, "--bam")
    program.add_arg(v, str, "-v")
    program.add_arg(vcf, str, "--vcf")
    program.add_arg(f, formats.Fasta, "-f")
    program.add_arg(fasta_reference, formats.Fasta, "--fasta-reference")
    program.add_arg(t, formats.Bed, "-t")
    program.add_arg(targets, formats.Bed, "--targets")
    program.add_arg(r, str, "-r")
    program.add_arg(region, str, "--region")
    program.add_arg(s, formats.TextFile, "-s")
    program.add_arg(samples, formats.TextFile, "--samples")
    program.add_arg(populations, formats.TextFile, "--populations")
    program.add_arg(A, formats.Bed, "-A")
    program.add_arg(cnv_map, formats.Bed, "--cnv-map")
    program.add_arg(trace, str, "--trace")
    program.add_arg(failed_alleles, formats.Bed, "--failed-alleles")
    program.add_arg(variant_input, formats.Vcf, "--variant-input")
    program.add_arg(l, bool, "-l")
    program.add_arg(only_use_input_alleles, bool, "--only-use-input-alleles")
    program.add_arg(haplotype_basis_alleles, formats.Vcf,
            "--haplotype-basis-alleles")
    program.add_arg(report_all_haplotype_alleles, bool,
            "--report-all-haplotype-alleles")
    program.add_arg(report_monorphic, bool, "--report-monorphic")
    program.add_arg(P, float, "-P")
    program.add_arg(pvar, float, "--pvar")
    program.add_arg(T, float, "-T")
    program.add_arg(theta, float, "--theta")
    program.add_arg(p, int, "-p")
    program.add_arg(ploidy, int, "--ploidy")
    program.add_arg(J, bool, "-J")
    program.add_arg(pooled_discrete, bool, "--pooled-discrete")
    program.add_arg(pooled_continuous, bool, "--pooled-continuous")
    program.add_arg(Z, bool, "-Z")
    program.add_arg(use_reference_allele, bool, "--use-reference-allele")
    program.add_args(reference_quality, int, 2, ",", "--reference-quality")
    program.add_arg(I, bool, "-I")
    program.add_arg(no_snps, bool, "--no-snps")
    program.add_arg(i, bool, "-i")
    program.add_arg(no_indels, bool, "--no-indels")
    program.add_arg(X, bool, "-X")
    program.add_arg(no_mnps, bool, "--no-mnps")
    program.add_arg(u, bool, "-u")
    program.add_arg(no_complex, bool, "--no-complex")
    program.add_arg(n, int, "-n")
    program.add_arg(use_best_n_alleles, int, "--use-best-n-alleles")
    program.add_arg(E, int, "-E")
    program.add_arg(max_complex_gap, int, "--max-complex-gap")
    program.add_arg(haplotype_length, int, "--haplotype-length")
    program.add_arg(min_repeat_length, int, "--min-repeat-length")
    program.add_arg(min_repeat_entropy, int, "--min-repeat-entropy")
    program.add_arg(no_partial_observations, bool, "--no-partial-observations")
    program.add_arg(O, bool, "-O")
    program.add_arg(dont_left_align_indels, bool, "--dont-left-align-indels")
    program.add_arg(_4, bool, "-4")
    program.add_arg(use_duplicate_reads, bool, "--use-duplicate-reads")
    program.add_arg(m, int, "-m")
    program.add_arg(min_mapping_quality, int, "--min-mapping-quality")
    program.add_arg(q, int, "-q")
    program.add_arg(min_base_quality, int, "--min-base-quality")
    program.add_arg(R, int, "-R")
    program.add_arg(min_supporting_allele_qsum, int,
            "--min-supportiong-allele-qsum")
    program.add_arg(Y, int, "-Y")
    program.add_arg(min_supporting_mapping_qsum, int,
            "--min-supportiong-mapping-qsum")
    program.add_arg(Q, int, "-Q")
    program.add_arg(mismatch_base_quality_threshold, int,
            "--mismatch-base-quality-threshold")
    program.add_arg(U, int, "-U")
    program.add_arg(read_mismatch_limit, int, "--read-mismatch-limit")
    program.add_arg(z, int, "-z")
    program.add_arg(read_max_mismatch_fraction, int,
            "--read-max-mismatch-fraction")
    program.add_arg(read_snp_limit, int, "--read-snp-limit")
    program.add_arg(e, int, "-e")
    program.add_arg(read_indel_limit, int, "--read-indel-limit")
    program.add_arg(_0, bool, "-0")
    program.add_arg(standard_filters, bool, "--standard-filters")
    program.add_arg(F, int, "-F")
    program.add_arg(min_alternate_fraction, int, "--min-alternate-fraction")
    program.add_arg(C, int, "-C")
    program.add_arg(min_alternate_count, int, "--min-alternate-count")
    program.add_arg(_3, int, "-3")
    program.add_arg(min_alternate_qsum, int, "--min-alternate-qsum")
    program.add_arg(G, int, "-G")
    program.add_arg(min_alternate_total, int, "--min-alternate-total")
    program.add_arg(min_coverage, int, "--min-coverage")
    program.add_arg(no_population_priors, bool, "--no-population-priors")
    program.add_arg(w, bool, "-w")
    program.add_arg(hwe_priors_off, bool, "--hwe-priors-off")
    program.add_arg(V, bool, "-V")
    program.add_arg(binomial_obs_priors_off, bool, "--binomial-obs-priors-off")
    program.add_arg(a, bool, "-a")
    program.add_arg(allele_balance_priors_off, bool,
            "--allele-balance-priors-off")
    program.add_arg(observation_bias, formats.TextFile, "--observation-bias")
    program.add_arg(base_quality_cap, int, "--base-quality-cap")
    program.add_arg(experimental_gls, bool, "--experimental-gls")
    program.add_arg(prob_contamination, float, "--prob-contamination")
    program.add_arg(contamination_estimates, formats.TextFile,
            "--contamination-estimates")
    program.add_arg(report_genotype_likelihood_max, bool,
            "--report-genotype-likelihood-max")
    program.add_arg(B, int, "-B")
    program.add_arg(genotyping_max_iterations, int,
            "--genotype-max-iterations")
    program.add_arg(genotyping_max_banddepth, int, "--genotype-max-banddepth")
    program.add_args(W, int, 2, ",", "-W")
    program.add_args(posterior_integration_limits, int, 2, ",",
            "--posterior-integration-limits")
    program.add_arg(exclude_unobserved_genotypes, bool,
            "--exclude-unobserved-genotypes")
    program.add_arg(S, int, "-S")
    program.add_arg(genotype_variant_threshold, int,
            "--genotype-variant-threshold")
    program.add_arg(j, bool, "-j")
    program.add_arg(use_mapping_quality, bool, "--use-mapping-quality")
    program.add_arg(H, bool, "-H")
    program.add_arg(harmonic_indel_quality, bool, "--harmonic-indel-quality")
    program.add_arg(D, int, "-D")
    program.add_arg(read_dependence_factor, int, "--read-dependence-factor")
    program.add_arg(genotype_qualities, bool, "--genotype-qualities")
    program.add_args(_in, formats.Bam, 1)
    return formats.Vcf(vcf, program)
    
