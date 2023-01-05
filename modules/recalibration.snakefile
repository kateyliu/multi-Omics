#MODULE: Indel and Base Realigner  by Sentieon
#import os
#from string import Template

_realigner_threads=32


def recalibration_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
    	ls.append("analysis/align/%s/%s.realigned.bam" % (sample,sample))
    	ls.append("analysis/align/%s/%s_prerecal_data.table" % (sample,sample))
    	ls.append("analysis/align/%s/%s_recalibrated.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s_postrecal_data.table" % (sample,sample))
        ls.append("analysis/align/%s/%s_recal.csv" % (sample,sample))
    if not config.get('tumor_only'): #Only run when we have normals
        for run in config['runs']:
            ls.append("analysis/corealignments/%s/%s_tn_corealigned.bam" % (run,run))
    return ls


rule recalibration_all:
    input:
        recalibration_targets
    benchmark: "benchmarks/recalibration/recalibration_all.txt"

rule Indel_realigner_sentieon:
    """indel realigner for uniquely  mapped reads"""
    input:
         bam="analysis/align/{sample}/{sample}.sorted.dedup.bam",
         bai="analysis/align/{sample}/{sample}.sorted.dedup.bam.bai",
    output:
         realignbam="analysis/align/{sample}/{sample}.realigned.bam",
         realignbai="analysis/align/{sample}/{sample}.realigned.bam.bai"
    message:
         "INDEL REALIGNER: indel realigner for  mapped reads"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        mills=config['Mills_indels'],
        g1000=config['G1000_indels'],
    group: "recalibration"
    threads: 8 #_realigner_threads
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Indel_realigner_sentieon.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.bam} --algo Realigner -k {params.mills} -k {params.g1000} {output.realignbam}"""

rule Base_recalibration_precal_sentieon:
    """base recalibration for realigned files"""
    input:
        realignbam="analysis/align/{sample}/{sample}.realigned.bam",
    output:
        prerecaltable="analysis/align/{sample}/{sample}_prerecal_data.table",
        recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam",
        recalibratedbai="analysis/align/{sample}/{sample}_recalibrated.bam.bai"
    message:
        " PRE BASE RECALIBRATION: base recalibration for  realigned files"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        mills= config['Mills_indels'],
        g1000= config['G1000_indels'],
    threads: 4 #_realigner_threads
    group: "recalibration"
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_precal_sentieon.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.realignbam} --algo QualCal -k {params.dbsnp} -k {params.mills} -k {params.g1000}  {output.prerecaltable} --algo ReadWriter {output.recalibratedbam}"""

rule Base_recalibration_postcal_sentieon:
    """post recalibration for realigned files"""
    input:
        recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam",
        recalibrated_bai="analysis/align/{sample}/{sample}_recalibrated.bam.bai",
        prerecaltable="analysis/align/{sample}/{sample}_prerecal_data.table"
    output:
        postrecaltable="analysis/align/{sample}/{sample}_postrecal_data.table"
    message:
        "POST BASE RECALIBRATION: post base recalibration for  realigned files"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        mills= config['Mills_indels'],
        g1000= config['G1000_indels'],
    threads: 16 #_realigner_threads
    group: "recalibration"
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_postcal_sentieon.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.recalibratedbam} -q {input.prerecaltable} --algo QualCal -k {params.dbsnp} -k {params.mills} -k {params.g1000}  {output.postrecaltable}"""

rule Base_recalibration_sentieon:
    """ recalibration for realigned files"""
    input:
        postrecaltable="analysis/align/{sample}/{sample}_postrecal_data.table",
        prerecaltable="analysis/align/{sample}/{sample}_prerecal_data.table"
    output:
        recalfile="analysis/align/{sample}/{sample}_recal.csv"
    message:
        "DIFF BASE RECALIBRATION: Difference in pre and post processing of realigned files"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
    threads: 1 #_realigner_threads
    group: "recalibration"
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_sentieon.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -t {threads} --algo QualCal --plot --before {input.prerecaltable} --after {input.postrecaltable} {output.recalfile}"""

def recal_corealignment_inputFn(wildcards):
    run = config['runs'][wildcards.run]
    normal = run[0]
    tumor = run[1]

    tmp = {}
    if not config.get('tumor_only'):
        tmp['normal'] = "analysis/align/%s/%s.realigned.bam" % (normal, normal)
        tmp['normal_recal'] = "analysis/align/%s/%s_prerecal_data.table" % (normal, normal)
    tmp['tumor'] = "analysis/align/%s/%s.realigned.bam" % (tumor, tumor)
    tmp['tumor_recal'] = "analysis/align/%s/%s_prerecal_data.table" % (tumor, tumor)
    return tmp

#NOTE: this rule shouldn't be called with tumor only runs, but we still have
# to define it.  It makes the params section ugly!
rule corealignment:
    input:
        unpack(recal_corealignment_inputFn)
    output:
        bam="analysis/corealignments/{run}/{run}_tn_corealigned.bam",
        bai="analysis/corealignments/{run}/{run}_tn_corealigned.bam.bai",
    params:
        #Change the files inputted depending on whether normal is available
        in_files = lambda wildcards,input: "-i %s" % input.tumor if config.get('tumor_only', False) else "-i %s -i %s" % (input.tumor, input.normal),
        recal_files = lambda wildcards,input: "-i %s" % input.tumor_recal if config.get('tumor_only', False) else "-q %s -q %s" % (input.tumor_recal, input.normal_recal),
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        #dbsnp= config['dbsnp'], #not used
        mills= config['Mills_indels'],
        g1000= config['G1000_indels'],
    threads: _realigner_threads
    group: "recalibration"
    benchmark:
        "benchmarks/recalibration/{run}/{run}.corealignment.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} {params.in_files} {params.recal_files} --algo Realigner -k {params.mills} -k {params.g1000} {output.bam}"""
