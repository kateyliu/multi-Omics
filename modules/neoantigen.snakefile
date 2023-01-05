#module: neoantigen prediction module using pvacseq pipeline (from pvactools)

_neoantigen_threads=64 #should be set to as max cores; w/ 64 runtime~=1hr

def neoantigen_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)
    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        tmp.append(sample_name)
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    return tmp


def neoantigen_getNormal(wildcards):
    return neoantigen_runsHelper(wildcards, 0)

def neoantigen_getTumor(wildcards):
    return neoantigen_runsHelper(wildcards, 1)

def neoantigen_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        tumor = config['runs'][run][1]
        print(tumor)
        ls.append("analysis/neoantigen/%s/MHC_Class_I/%s.filtered.tsv" % (run,tumor))
    return ls


def neoantigen_get_output_keys(wildcards):
    #ls = ['classI_ranked_epitopes','classI_all_epitopes','classII_ranked_epitopes','classII_all_epitopes',
       #               'combined_ranked_epitopes','combined_all_epitopes']
    ls = ['classI_ranked_epitopes','classI_all_epitopes','combined_all_epitopes']
    #if 'neoantigen_run_classII' in config and config['neoantigen_run_classII']:
        #ls.extend(['classII_ranked_epitopes','classII_all_epitopes',
                   #'combined_ranked_epitopes','combined_all_epitopes'])
    return ls

def neoantigen_getNeoantigenList_helper(wildcards):
    ls = []
    run = wildcards.run
    tumor = neoantigen_getTumor(wildcards)[0]
    ls.append("analysis/neoantigen/%s/combined/%s.filtered.tsv" % (run,tumor))
    #ls.append("analysis/somatic/%s/%s_tnscope.filter.vep.vcf.gz" % (run,tumor))
    return ls

def getPvacseqOut(wildcards):
    """returns tuple (run, tumor sample name)"""
    run = wildcards.run
    tumor = config['runs'][run][1]
    return(run,tumor)

def getVCF_file(wildcards):
    """IF there is expression data available, return the expression-added VCF
    otherwise return the neoantigen prepared vcf"""

    run = wildcards.run
    tumor_sample = config['runs'][run][1]
    ret = "analysis/somatic/%s/%s_tnscope.filter.vep.vcf" % (run,run)
    return ret

def getTumorHLA(wildcards):
    """get the optitype results file for the tumor sample"""
    run = wildcards.run
    tumor = config['runs'][run][1]
    ls = ["analysis/optitype/%s/%s_result.tsv" % (tumor, tumor)]
    #if 'neoantigen_run_classII' in config and config['neoantigen_run_classII']:
        #ls.append("analysis/hlahd/%s/result/%s_final.result.txt" % (tumor,tumor))
    return ls

def parseHLA(hla_files):
    """Given an optitypes '_results.tsv' file; parses the HLA A, B, C
    and returns these as a comma-separated string (for pvacseq) input

    NOTE: cureently the optitype results.tsv looks somthing like this:
    	A1	A2	B1	B2	C1	C2	Reads	Objective
    0					C*06:04	C*06:04	4.0	3.99
    **So were' going to parse cols 1-6 and return that"""
    #Internal helper fn
    def _cleanAllele(s):
        """takes a strine HLA-DRB1*13:02:01 to DRB1*13:02"""
        #Remove the HLA-
        s = s.split("-")[1]
        #Remove the last :01
        s = ":".join(s.split(":")[:-1])
        #print(s)
        return s

    #CATCH when the HLA does not exist yet
    #print(optitype_out_file)
    optitype_out_file = hla_files[0]
    if not os.path.exists(optitype_out_file):
        #print("WES WARNING: %s is not found!" % optitype_out_file)
        return ""

    f = open(optitype_out_file)
    hdr = f.readline().strip().split("\t") #ignore for now
    classI = f.readline().strip().split("\t")[1:7] #want first 6 cols
    #FOR classI alleles, prepend a HLA to each of them
    classI = ["HLA-%s" % a for a in classI if a]
    #print(classI)
    f.close()

    #check for xhla file
#    classII = []
#    if 'neoantigen_run_classII' in config and config['neoantigen_run_classII'] and len(hla_files) > 1:
#        classII_out_file = hla_files[1]

        #SUB-in this section with a call to report_neoantigens_hla.parseHLA_HD
        #PARSE hlahd txt file...
 #       c2_alleles = {}
        #NOTE read only A,B,C,DRB1,DQA1,DQB1,DPA1,DPB1 (first 8 lines)
        #REST: DMA,DMB,DOA,DOB,DRA,DRB2,DRB3,DRB4,DRB5,DRB6,DRB7,DRB8,DRB9,
        #DPA2,E,???,G,H,J,K,L,???,V,???,Y
#        if os.path.exists(classII_out_file):
#            f = open(classII_out_file)
#            for i in range(8): #first 8 lies
#                l = f.readline().strip()
#                if i > 2 and not l.startswith('Could'):#ksip 'Couldn't read result file.' and A,B,C alleles
#                    tmp = l.split("\t")
#                    #Coerce HLA-DRB1*13:02:01 to DRB1*13:02
#                    c2_alleles[tmp[0]] = [_cleanAllele(a) for a in tmp[1:3] if a != "-"]
#            f.close()
#        #END SUB
#            #print(c2_alleles)
#            #COMPOSE the classII alleles- DRP1 just add them
#            for a in c2_alleles.get('DRB1', []):
#                classII.append(a)
#            for a in c2_alleles.get('DQA1', []): #DQ
#                #match a with each DQB1
#                for b in c2_alleles.get('DQB1', []):
#                    classII.append("%s-%s" % (a,b))
#            for a in c2_alleles.get('DPA1', []): #REPEAT for DP
#                for b in c2_alleles.get('DPB1', []):
#                    classII.append("%s-%s" % (a,b))
            #print(classII)

#    if classII:
#        classI.extend(classII)
#    #NOTE: NOW classI has all hla alleles (including classII if opted for)
    hla = ",".join(["%s" % a for a in classI if a])
    print(hla)
    return hla




rule neoantigen_all:
    input:
        neoantigen_targets
    benchmark: "benchmarks/neoantigen/neoantigen_all.txt"


rule neoantigen_pvacseq:
        """NOTE: neoantigen's pvacseq is not available on CONDA
        MUST either be installed in base system/docker container"""
        input:
            vcf=getVCF_file,
            hla=getTumorHLA,
        output:
            filtered="analysis/neoantigen/{run}/MHC_Class_I/{tumor}.filtered.tsv",
            all_epitopes="analysis/neoantigen/{run}/MHC_Class_I/{tumor}.all_epitopes.tsv",
        params:
            normal = lambda wildcards: "--normal-sample-name %s" % config['runs'][wildcards.run][0] if not config.get('tumor_only') else "",
            tumor = lambda wildcards: config['runs'][wildcards.run][1],
            iedb = config['neoantigen_iedb'],
            HLA = lambda wildcards,input: parseHLA(input.hla),
            callers=config.get('neoantigen_callers','NetMHCpan'),
            epitope_lengths_cls1=config.get('neoantigen_epitope_lengths_cls1', '8,9,10,11'),
            epitope_lengths_cls2=config.get('neoantigen_epitope_lengths_cls2', '12,13,14,15,16,17,18'),
            output_dir = lambda wildcards: "%sanalysis/neoantigen/%s/" % (config['remote_path'], wildcards.run),
        threads: 64 #_neoantigen_threads
        group: "neoantigen"
        log: "analysis/logs/neoantigen/{run}/{tumor}.neoantigen_pvacseq.log"
        benchmark:
            "benchmarks/neoantigen/{run}/{tumor}.neoantigen_pvacseq.txt"
        shell:
            """pvacseq run {input.vcf} {params.tumor} {params.HLA} {params.callers} {params.output_dir} -e1 {params.epitope_lengths_cls1} -t {threads} {params.normal} --iedb-install-directory {params.iedb} 2> {log}"""


