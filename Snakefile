################################################################################
## CNV Pipeline
## Copy Number Variant (CNV) discovery with Sentieon CNV
##
## Authors: Massimiliano Volpe, Jyotirmoy Das
## Email: massimiliano.volpe@liu.se, jyotirmoy.das@liu.se
## Date: 11/02/2021
## Developed on behalf of the Bioinformatics Core Facility, Linköping University
##
## Rules:

################################################################################


# Functions -------------------------------------------------------------------
def id_maker(path, sep, sample, suffix, d):
    f = "".join([path, sep, sample, suffix])
    l = subprocess.check_output("zcat " + f + " | head -n 1 | cut -d ':' -f 3,4 | sed  s/:/./g | sed 's/@//'", shell=True).strip().decode()
    d[sample] = l
    return(d)

def pu_maker(path, sep, sample, suffix, d):
    f = "".join([path, sep, sample, suffix])
    l = subprocess.check_output("zcat " + f + " | head -n 1 | cut -d ':' -f 3,4,10 | sed  s/:/./g | sed 's/@//'", shell=True).strip().decode()
    d[sample] = l
    return(d)


# Globals ---------------------------------------------------------------------
import subprocess

configfile:
    "config.json"

#workdir:
#    config['workdir']


R1SUFFIX = config['R1_suffix']
R2SUFFIX = config['R2_suffix']

SAMPLES, = glob_wildcards(config['normal'] + "/{sample}" + R1SUFFIX)
RESULTS = config['workdir'] + '/{sample}/'
BAMS = RESULTS + 'bams/'
LOGS = RESULTS + 'logs/'
METRICS = RESULTS + 'metrics/'
PLOTS = RESULTS + 'plots/'
MARKDUP = RESULTS + 'markdup/'
RECAL = RESULTS + 'baserecal/'

ref = config['reference']
dbsnp = config['dbsnp']
known_Mills_indels = config['known_Mills_indels']
known_1000G_indels = config['known_1000G_indels']
interval = config['interval']

id = {}
pu = {}

for sample in SAMPLES:
    id_maker(config['normal'], '/', sample, R1SUFFIX, id)
    pu_maker(config['normal'], '/', sample, R1SUFFIX, pu)


# Rules -----------------------------------------------------------------------
rule all:
    input:
        #normal_bam = expand(BAMS + "{sample}.normal_sorted.bam", sample=SAMPLES),
        #normal_ml = expand(LOGS + '{sample}.normal_metrics.log', sample=SAMPLES),
        #normal_deduped = expand(BAMS + '{sample}.normal_deduped.bam', sample=SAMPLES),
        #normal_rdt = expand(RECAL + '{sample}.normal_recal_data.table', sample=SAMPLES),
        bam = expand(BAMS + '{sample}.normal_recaled.bam', sample=SAMPLES),
        cov = config['workdir'] + '/coverage_file',
        pon = config['workdir'] + '/PoN.log'


rule mapping_normal:
    input:
        R1 = config['normal'] + "/{sample}" + R1SUFFIX,
        R2 = config['normal'] + "/{sample}" + R2SUFFIX,
        ref = config['reference']
    output:
        sam = temp(BAMS + '{sample}.normal.sam'),
        bam = temp(BAMS + '{sample}.normal_sorted.bam')
    log:
        bwa = LOGS + '{sample}.normal_bwa.log',
        sort = LOGS + '{sample}.normal_sort.log'
    params:
        #ID = subprocess.check_output("zcat {input.R1} | head -n 1 | cut -d ':' -f 1,4 | sed  s/:/./ | sed 's/@//'", shell=True).strip().decode(),
        #R = "@RG\\tID:" + config["normal_group"] + "\\tSM:{sample}_normal\\tPL:" + config["platform"]
        K = 10000000,
        ID = lambda wildcards: id[wildcards.sample],
        PU = lambda wildcards: pu[wildcards.sample],
        SM = "{sample}_normal",
        PL = config["platform"]
    threads:
        48 # set the maximum number of available cores
    shell:
        # sentieon bwa mem -M -R '{params.R}' -t {threads} -K {params.K} -o {output.sam} {input.ref} {input.R1} {input.R2} >> {log.bwa} 2>&1
        """
        sentieon bwa mem -M -R '@RG\\tID:{params.PU}\\tSM:{params.SM}\\tPL:{params.PL}' -t {threads} -K {params.K} -o {output.sam} {input.ref} {input.R1} {input.R2} >> {log.bwa} 2>&1
        sentieon util sort -r {input.ref} -i {output.sam} -o {output.bam} -t {threads} --sam2bam >> {log.sort} 2>&1
        """


rule metrics_normal:
    input:
        bam = rules.mapping_normal.output.bam,
        ref = config['reference']
    output:
        mqm = METRICS + '{sample}.normal_mq_metrics.txt',
        qdm = METRICS + '{sample}.normal_qd_metrics.txt',
        gcs = METRICS + '{sample}.normal_gc_summary.txt',
        gcm = METRICS + '{sample}.normal_gc_metrics.txt',
        aln = METRICS + '{sample}.normal_aln_metrics.txt',
        ism = METRICS + '{sample}.normal_is_metrics.txt',
        gcp = PLOTS + '{sample}.normal_gc-report.pdf',
        qdp = PLOTS + '{sample}.normal_qd-report.pdf',
        mqp = PLOTS + '{sample}.normal_mq-report.pdf',
        isp = PLOTS + '{sample}.normal_is-report.pdf'
    log:
        LOGS + '{sample}.normal_metrics.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        sentieon driver -r {input.ref} -t {threads} -i {input.bam} --algo MeanQualityByCycle {output.mqm} --algo QualDistribution {output.qdm} --algo GCBias --summary {output.gcs} {output.gcm} --algo AlignmentStat --adapter_seq '' {output.aln} --algo InsertSizeMetricAlgo {output.ism} >> {log} 2>&1
        sentieon plot GCBias -o {output.gcp} {output.gcm}
        sentieon plot QualDistribution -o {output.qdp} {output.qdm}
        sentieon plot MeanQualityByCycle -o {output.mqp} {output.mqm}
        sentieon plot InsertSizeMetricAlgo -o {output.isp} {output.ism}
        """


rule markdup_normal:
    input:
        bam = rules.mapping_normal.output.bam,
        ref = config['reference']
    output:
        ns = MARKDUP + '{sample}.normal_score.txt',
        dm = MARKDUP + '{sample}.normal_dedup_metrics.txt',
        bam = temp(BAMS + '{sample}.normal_deduped.bam'),
        cm = MARKDUP + '{sample}.normal_coverage_metrics'
    log:
        LOGS + '{sample}.normal_dedup.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        sentieon driver -t {threads} -i {input.bam} --algo LocusCollector --fun score_info {output.ns} >> {log} 2>&1
        sentieon driver -t {threads} -i {input.bam} --algo Dedup --rmdup --score_info {output.ns} --metrics {output.dm} {output.bam} >> {log} 2>&1
        sentieon driver -r {input.ref} -t {threads} -i {output.bam} --algo CoverageMetrics {output.cm} >> {log} 2>&1
        """


rule baserecal_normal:
    input:
        bam = rules.markdup_normal.output.bam,
        ref = config['reference']
    output:
        rdt = RECAL + '{sample}.normal_recal_data.table',
        post = RECAL + '{sample}.normal_recal_data.table.post',
        recal = RECAL + '{sample}.normal_recal.csv',
        rp = PLOTS + '{sample}.normal_recal_plots.pdf',
        bam = BAMS + '{sample}.normal_recaled.bam'
    log:
        LOGS + '{sample}.normal_recal.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        sentieon driver -r {input.ref} -t {threads} -i {input.bam} --algo QualCal -k {dbsnp} -k {known_Mills_indels} -k {known_1000G_indels} {output.rdt} >> {log} 2>&1
        sentieon driver -r {input.ref} -t {threads} -i {input.bam} -q {output.rdt} --algo QualCal -k {dbsnp} -k {known_Mills_indels} -k {known_1000G_indels} {output.post} --algo ReadWriter {output.bam} >> {log} 2>&1
        sentieon driver -t {threads} --algo QualCal --plot --before {output.rdt} --after {output.post} {output.recal} >> {log} 2>&1
        sentieon plot QualCal -o {output.rp} {output.recal}
        """


l = []
for i in expand(BAMS + '{sample}.normal_recaled.bam', sample=SAMPLES): l.append('-i ' + i)
#print(l)


rule coverage_file:
    input:
        bam = expand(BAMS + '{sample}.normal_recaled.bam', sample=SAMPLES),
        ref = config['reference']
    output:
        cov = config['workdir'] + '/coverage_file',
    log:
        config['workdir'] + '/coverage_file.log'
    params:
        i = l
    threads:
        48 # set the maximum number of available cores
    shell:
        # sentieon driver -t {threads} -r {input.ref} {input.bam} --algo CNV --target BED_FILE [--target_padding PADDING  --coverage COVERAGE_FILE] --create_coverage {output.cov}
        """
        sentieon driver -t {threads} -r {input.ref} {params.i} --algo CNV --target {interval} --create_coverage {output.cov} >> {log} 2>&1
        """


rule create_PoN:
    input:
        bam = expand(BAMS + '{sample}.normal_recaled.bam', sample=SAMPLES),
        ref = config['reference']
    output:
        pon = config['workdir'] + '/PoN',
    log:
        config['workdir'] + '/PoN.log'
    params:
        i = l
    threads:
        48 # set the maximum number of available cores
    shell:
        # sentieon driver -t {threads} -r {input.ref} {input.bam} --algo CNV --target BED_FILE [--target_padding PADDING  --coverage COVERAGE_FILE] --create_coverage {output.cov}
        """
        sentieon driver -t {threads} -r {input.ref} {params.i} --algo CNV --target {interval} --create_pon {output.pon} >> {log} 2>&1
        """


############
# COMMENTS #
############
# flowcell:lane is actually problematic since this error could arise:
# Readgroup H3M2WBGXY.1 with different attributes is present in multiple bam files: /media/Data2/mouna_testing_data/results/PoN/342_S8/bams/342_S8.normal_recaled.bam, 
# /media/Data2/mouna_testing_data/results/PoN/339_S7/bams/339_S7.normal_recaled.bam, /media/Data2/mouna_testing_data/results/PoN/265_S2/bams/265_S2.normal_recaled.bam, 
# /media/Data2/mouna_testing_data/results/PoN/293_S3/bams/293_S3.normal_recaled.bam, /media/Data2/mouna_testing_data/results/PoN/343_S9/bams/343_S9.normal_recaled.bam, 
# /media/Data2/mouna_testing_data/results/PoN/318_S5/bams/318_S5.normal_recaled.bam, /media/Data2/mouna_testing_data/results/PoN/311_S4/bams/311_S4.normal_recaled.bam, 
# /media/Data2/mouna_testing_data/results/PoN/337_S6/bams/337_S6.normal_recaled.bam, /media/Data2/mouna_testing_data/results/PoN/256_S1/bams/256_S1.normal_recaled.bam.
# Error: Invalid input BAM files
# fixerd using PU as ID
