#module: somatic tumor calling using TNscope


def somatic_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(rr)
    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        tmp.append(sample_name)
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    return tmp


def somatic_getNormal(wildcards):
    return somatic_runsHelper(wildcards, 0)

def somatic_getTumor(wildcards):
    return somatic_runsHelper(wildcards, 1)


def somatic_targets(wildcards):
    ls = []
    #center = config.get('cimac_center', 'broad') #Try to get center, default broad
    for run in config['runs']:
        ls.append("analysis/somatic/%s/%s_tnscope.output.vcf.gz" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.filter.vcf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.filter.vep.vcf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.filter.vcf.gz" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.filter.vep.vcf.gz" % (run,run))
        #ls.append("analysis/somatic/%s/%s_tnscope.filter.vcf.gz.tbi" % (run,run))
        #ls.append("analysis/somatic/%s/%s_tnscope.filter.vep.vcf.gz.tbi" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.filter.vep.vcf.maf" % (run,run))
    return ls

rule somatic_all:
    input:
        somatic_targets
    benchmark: "benchmarks/somatic/somatic_all.txt"


rule somatic_calling_tumor_TNscope:
    input:
      corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
    output:
      tnscopevcf="analysis/somatic/{run}/{run}_tnscope.output.vcf.gz",
	  #tnscopevcf_tbi="analysis/somatic/{run}/{run}_tnscope.output.vcf.gz.tbi"
    params:
      index=config['genome_fasta'],
      sentieon_path=config['sentieon_path'],
      dbsnp= config['dbsnp'],
      tumor = lambda wildcards: config['runs'][wildcards.run][1],
      normal = lambda wildcards: config['runs'][wildcards.run][0],
      trim_soft_clip = "--trim_soft_clip" if config.get("trim_soft_clip", False) else "",
    threads: 18 #_somatic_threads
    priority: 50
    group: "somatic"
    benchmark:
      "benchmarks/somatic/{run}/{run}.somatic_calling_TNscope.txt"
    shell:
      """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads}  -i {input.corealignedbam} --algo TNscope --tumor_sample {params.tumor} --normal_sample {params.normal} --dbsnp {params.dbsnp} {params.trim_soft_clip} {output.tnscopevcf}"""

rule filter_raw_vcf:
    """General rule to filter the three different types of vcf.gz files"""
    input:
      "analysis/somatic/{run}/{run}_tnscope.output.vcf.gz"
    output:
      "analysis/somatic/{run}/{run}_tnscope.filter.vcf"
    params:
      index=config['genome_fasta'],
      sentieon_path=config['sentieon_path'],
      tumor=lambda wildcards: somatic_getTumor(wildcards),
      normal= lambda wildcards: somatic_getNormal(wildcards),
      vcf_bin_path="%s/bin/" % config['vcf_root'],
    group: "somatic"
    benchmark:
      "benchmarks/somatic/{run}/{run}.tnscope_filter_raw_vcf.txt"
    shell:
      """{params.vcf_bin_path}vcftools --gzvcf {input} --remove-filtered-all --recode --stdout > {output}"""
            #"""vcftools --gzvcf {input} --remove-filtered-all --recode --stdout > {output}"""

rule vcfVEP:
    """Rule to annotate vcf files with vep"""
    input:
      "analysis/somatic/{run}/{run}_tnscope.filter.vcf"
    output:
      "analysis/somatic/{run}/{run}_tnscope.filter.vep.vcf"
    params:
      vep_data=config['vep_data'],
      #vep_synonyms=config['vep_synonyms'],
      vep_plugins=config['vep_plugins'],
      cache_version='86',
      gdc_fasta=config['genome_fasta'],
      vcf_bin_path="%s/bin/" % config['vcf_root'],
      path="set +eu;source activate %s" % config['vep_root']
	#path= "set +eu;source activate /liulab/jinwang/pipeline_test/rnaseq_pipeline/envs/vep_env" #added jins env
    benchmark:
      "benchmarks/somatic/{run}/{run}.tnscope_vcfVEP.txt"
    group: "somatic"
    shell:
      #"{params.vcf_bin_path}vep --i {input} --dir_cache={params.vep_data} --synonyms {params.vep_synonyms} --vcf -o {output} --offline --hgvs --fa {params.gdc_fasta} --format vcf"
            #"vep --i {input} --dir_cache={params.vep_data} --synonyms {params.vep_synonyms} --vcf -o {output} --offline --hgvs --fa {params.gdc_fasta} --format vcf"
       #from jacob--"""{params.vcf_bin_path}vep --input_file {input} --output_file {output} --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta {params.gdc_fasta} --offline --cache --dir_cache {params.vep_data} --cache_version 84 --plugin Downstream --plugin Wildtype --dir_plugins {params.vep_plugins} --pick"""
        """{params.path};vep --input_file {input} --output_file {output} --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta {params.gdc_fasta} --offline --cache --dir_cache {params.vep_data} --cache_version 84 --plugin Downstream --plugin Wildtype --plugin Frameshift  --dir_plugins {params.vep_plugins} --pick"""

rule somatic_tabix_filtered_vcf_gz:
    """Prepping the files filtered.vcf file for somatic_getExonic_mutations
    bgzip-ing and tabix"""
    input:
        filtered_in="analysis/somatic/{run}/{run}_tnscope.filter.vcf",
        vep_in="analysis/somatic/{run}/{run}_tnscope.filter.vep.vcf"
    output:
        filtered_vcf="analysis/somatic/{run}/{run}_tnscope.filter.vcf.gz",
        vep_vcf="analysis/somatic/{run}/{run}_tnscope.filter.vep.vcf.gz",
        filtered_out="analysis/somatic/{run}/{run}_tnscope.filter.vcf.gz.tbi",
        vep_out="analysis/somatic/{run}/{run}_tnscope.filter.vep.vcf.gz.tbi"
    benchmark:
        "benchmarks/somatic/{run}/{run}_tnscope.tabix_filtered_vcf_gz.txt"
    group: "somatic"
    #conda: "../envs/somatic_vcftools.yml"
    shell:
        "bgzip -c {input.filtered_in} > {output.filtered_vcf}"
        "&& bgzip -c {input.vep_in} > {output.vep_vcf}"
        "&& tabix -p vcf {output.filtered_vcf}"
        "&& tabix -p vcf {output.vep_vcf}"

rule vcf2maf:
    """General rule to convert the different vcf files into maf"""
    input:
      vep="analysis/somatic/{run}/{run}_tnscope.filter.vep.vcf",
    output:
      "analysis/somatic/{run}/{run}_tnscope.filter.vep.vcf.maf"
    params:
      vep_index=config['vep_fasta'],
     # vep_custom_enst= config['vep_custom_enst'],
      vep_assembly=config['vep_assembly'],
      #vep_filter= config['vep_filter'],
      buffer_size=config['vcf2maf_bufferSize'],
      vcf_bin_path="%s/bin/" % config['vcf_root'],
      tumor= lambda wildcards: somatic_getTumor(wildcards),
      normal= lambda wildcards: somatic_getNormal(wildcards),
    benchmark:
      "benchmarks/somatic/{run}/{run}.tnscope_vcf2maf.txt"
    log:
      "analysis/logs/somatic/{run}/{run}.tnscope_vcf2maf.log.txt"
    group: "somatic"
    conda: "../envs/vcf.yml"
    shell:
        """{params.vcf_bin_path}vcf2maf.pl --input-vcf {input.vep} --output-maf {output} --ref-fasta {params.vep_index} --tumor-id {params.tumor} --normal-id {params.normal} --ncbi-build {params.vep_assembly} --buffer-size {params.buffer_size} --inhibit-vep 1 2> {log}"""
