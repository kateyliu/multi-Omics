#!/usr/bin/env python

from genson import SchemaBuilder

###############################################################################
# Constants
###############################################################################
_int = 777
_float = 7.0
_string = "TGBTG"
_int_arr = [_int,_int,_int]
_file_path = "/the/only/way"

###############################################################################
# SAMPLE level definition
###############################################################################
alignment = {'total_reads': _int,
             'mapped_reads': _int,
             'dedup_reads': _int,
             'gc_content': _int_arr, #NOTE: could be floats
             'quality_score': _int_arr} #NOTE: could be floats

coverage = {'total_reads': _int,
            'mean_depth': _float,
            "q1_depth": _float,
            "median_depth": _float,
            "q3_depth": _float,
            "percent_bases_gt_50": _float}

hla = { "A-1": _string,
        "A-2": _string,
        "B-1": _string,
        "B-2": _string,
        "C-1": _string,
        "C-2": _string,
        "DPB1-1": _string, #OPTIONAL
        "DPB1-2": _string, #OPTIONAL
        "DQB1-1": _string, #OPTIONAL
        "DQB1-2": _string, #OPTIONAL
        "DRB1-1": _string, #OPTIONAL
        "DRB1-2": _string} #OPTIONAL

sample = {'id': _string,
          'raw_file1': _file_path, #bucket path, either fastq or bam
          'raw_file2': _file_path, #bucket path, in case of PE
          'dedup_bam_file': _file_path,
          'alignment': alignment,
          'coverage': coverage,
          'hla': hla,
          'msi_sensor': _float}

###############################################################################
# END SAMPLE level definition
###############################################################################

###############################################################################
# SOMATIC level definition
###############################################################################

mutation_results = {"total": _int,
                    "snp": _int,
                    "insertion": _int,
                    "deletion": _int}

transition_row = {"A": _int,
                  "C": _int,
                  "G": _int,
                  "T": _int}
tmb = {'tumor': _int, 'normal': _int, "common": _int, "overlap": _float}
func_summary = {'missense': mutation_results, 'nonsense': mutation_results,
                'silent': mutation_results}

somatic_results = {'filtered_vcf_file': _file_path,
                   'filtered_maf_file': _file_path,
                   'tmb': tmb,
                   'mutation_summary': mutation_results,
                   'functional_summary': func_summary,
                   'trinucleotide_matrix': _int_arr,
                   'transition_matrix': {'A': transition_row,
                                         'C': transition_row,
                                         'G': transition_row,
                                         'T': transition_row}}
###############################################################################
# END SOMATIC level definition
###############################################################################

###############################################################################
# RUN level definition
###############################################################################

neoantigen_row = {"Gene_Name": _string,
                  "Mutation_Protein": _string,
                  "Position": _int,
                  "HGVSc": _string,
                  "HGVSp": _string,
                  "HLA_Allele": _string,
                  "MT_Epitope_Seq": _string,
                  "MT_IC50": _float,
                  "WT_IC50": _float,
                  "Fold_Change": _float,
                  "Tumor_DNA_Depth": _int,
                  "Tumor_DNA_VAF": _float,
                  "Score": _float}

copy_number = {'clonality': _float,
               'purity': _float,
               'ploidy': _float,
               'dipLogR': _float,
               'cnv_file': _file_path, #path to sequenza result txt file
               'cnv_plot_file': _file_path}

run = {'id': _string,
       'tumor': sample,
       'normal': sample,
       'copy_number': copy_number,
       'somatic': somatic_results,
       'neoantigen': [neoantigen_row, neoantigen_row, neoantigen_row],
       'neoantigen_file': _file_path, #path to pvacseq filtered.condensed.ranked.tsv file
       }

###############################################################################
# END RUN level definition
###############################################################################

builder = SchemaBuilder(schema_uri="http://json-schema.org/draft-07/schema#")
builder.add_object(run)
print(builder.to_json(indent=3))
