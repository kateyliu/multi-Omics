# module:Precision HLA typing from next-generation sequencing data by Optitype and Polysolver
_optitype_threads=64

def hla_targets(wildcards):
    """ Generates the targets for this module"""
    ls = []
    for sample in config['samples']:
        ls.append("analysis/optitype/%s/%s_result.tsv" % (sample,sample))
        ls.append("analysis/optitype/%s/%s_coverage_plot.pdf" % (sample,sample))
        ls.append("analysis/optitype/%s/%s.sorted.chr6.bam" % (sample,sample))
        ls.append("analysis/optitype/%s/%s.sorted.chr6.end1.fastq" % (sample,sample))
        ls.append("analysis/optitype/%s/%s.sorted.chr6.end2.fastq" % (sample,sample))
    #    ls.append("analysis/hlahd/%s/result/%s_final.result.txt" % (sample,sample))
    return ls

rule hla_all:
    input:
        hla_targets
    benchmark: "benchmarks/optitype/hla_all.txt"

rule optitype_extract_chr6:
    """Extract chr6 by sambamba"""
    input:
        in_sortbamfile = "analysis/align/{sample}/{sample}.sorted.dedup.bam"
    output:
        chr6sortbamfile = "analysis/optitype/{sample}/{sample}.sorted.chr6.bam"
    threads: 16 #_optitype_threads
    group: "optitype"
    conda: "../envs/optitype.yml"
    benchmark:
        "benchmarks/optitype/{sample}/{sample}.optitype_extract_chr6.txt"
    shell:
        """sambamba view -t {threads} -f bam -h {input.in_sortbamfile} chr6 > {output.chr6sortbamfile}"""

rule optitype_index_chr6bam:
    """index chr6bam"""
    input:
        "analysis/optitype/{sample}/{sample}.sorted.chr6.bam"
    output:
        "analysis/optitype/{sample}/{sample}.sorted.chr6.bam.bai"
    threads: 1 #_optitype_threads
    group: "optitype"
    conda: "../envs/optitype.yml"
    benchmark:
        "benchmarks/optitype/{sample}/{sample}.optitype_index_chr6bam.txt"
    shell:
        """sambamba index -t {threads} {input}"""

rule optitype_bamtofastq:
    """Convert the sorted.chr6.bam file to fastq by samtools"""
    input:
        in_sortchr6bamfile = "analysis/optitype/{sample}/{sample}.sorted.chr6.bam",
        in_index = "analysis/optitype/{sample}/{sample}.sorted.chr6.bam.bai"
    output:
        chr6fastqfile1 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end1.fastq",
        chr6fastqfile2 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end2.fastq"
    group: "optitype"
    conda: "../envs/optitype.yml"
    log: "analysis/logs/optitype/{sample}/{sample}.optitype_bamtofastq.log"
    benchmark:
        "benchmarks/optitype/{sample}/{sample}.optitype_bamtofastq.txt"
    shell:
        """samtools fastq -@ 2 -1 {output.chr6fastqfile1} -2 {output.chr6fastqfile2} {input.in_sortchr6bamfile} 2> {log}"""

rule optitype_hlatyping:
    """Precision HLA typing from next-generation sequencing data by
    OptiType This will produce a time-stamped directory inside the
    specified outputn directory containing a CSV file with the predicted
    optimal (and if enumerated, sub-optimal)HLA genotype, and a pdf file
    containing a coverage plot of the predicted alleles for diagnostic
    purposes"""
    input:
        in_chr6fastqfile1 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end1.fastq",
        in_chr6fastqfile2 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end2.fastq"
    params:
        #PathtoOptiType = config['conda_path'],
        name=lambda wildcards: wildcards.sample,
        output_dir=lambda wildcards: "%sanalysis/optitype/%s/" % (config['remote_path'], wildcards.sample),
        #outputname = lambda wildcards: wildcards.sample
        path="set +u; source activate %s" % config['optitype_root'],
        optitype_config="wes/static/optitype/config.ini",
    output:
        HLAgenotype = "analysis/optitype/{sample}/{sample}_result.tsv",
        Coverageplot = "analysis/optitype/{sample}/{sample}_coverage_plot.pdf"
    group: "optitype"
    conda: "../envs/optitype.yml"
    benchmark:
        "benchmarks/optitype/{sample}/{sample}.optitype_hlatyping.txt"
    shell:
        """{params.path}; OptiTypePipeline.py -i {input.in_chr6fastqfile1} {input.in_chr6fastqfile2} --dna -v -o {params.output_dir} -p {params.name} --config {params.optitype_config}"""


#rule hlahd:
#    """calculate hlatyping by hla-hd"""
#    input:
#        chr6fastqfile1 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end1.fastq",
#        chr6fastqfile2 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end2.fastq"
#    output:
#        "analysis/hlahd/{sample}/result/{sample}_final.result.txt"
#    threads: _optitype_threads
#    group: "hlahd"
#    params:
#        name=lambda wildcards: wildcards.sample,
#        output_dir=lambda wildcards: "%sanalysis/hlahd/" % config['remote_path'],
#        min_lengh = 50, #hla-hd param -m which sets the min read length- fixed b/c we don't expect any libraries with shorter than 50bp
#        cut_perc = 0.95, #hla-hd param -c, if mistmatch how much to cut back
#        freq_data = config['hlahd_freq_data'],
#        split_file = config['hlahd_split_file'],
#        dictionary = config['hlahd_dictionary'],
#    log: "analysis/logs/hlahd/{sample}/{sample}.hlahd.log"
#    benchmark:
#        "benchmarks/hlahd/{sample}/{sample}.hlahd.txt"
#    shell:
#        """hlahd.sh -m {params.min_lengh} -c {params.cut_perc} -t {threads} -f {params.freq_data} {input.chr6fastqfile1} {input.chr6fastqfile2} {params.split_file} {params.dictionary} {params.name} {params.output_dir} 2> {log}"""
