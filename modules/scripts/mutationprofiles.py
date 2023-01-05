""" Visualization of WES harmonization 2
Usage: cidc-vs /path/to/files /path/to/output

Input:

    Path of folder for all files

    Path of folder for output figures

Files nomenclature:

    Center_source_sampleID_tools_feature
    Options:
        Center: Broad, MDACC
        tools: Sentieon, GISTIC, Pyclone, Facets
        source: FFPE, Frozen
        feature: Mutation.Germline, Mutation.Somatic, Clonality, Purity, CNV

"""
from .preprocess import *
from .definition import *
from .plot import *
import warnings
from matplotlib.backends.backend_pdf import PdfPages

BASE_DIR = os.path.dirname(os.path.realpath(__file__))+'/'
REF = os.path.join(BASE_DIR, 'data/REF')
class Compare(object):
    __slots__ = [
                "somatic_mutation_maf_dir", "somatic_mutation_vcf_dir",
                 "ref_fasta", "cancer",
                ]

    def __init__(out_dir,
                somatic_mutation_maf_dir,somatic_mutation_vcf_dir, ref_fasta,cancer):
        self.cancer = cancer
        self.ref_fasta = ref_fasta
        self.out_dir = out_dir
        self.somatic_mutation_maf_dir = somatic_mutation_maf_dir
        self.somatic_mutation_vcf_dir = somatic_mutation_vcf_dir
        self.DPI = 100
        self.fmt='png'

        def somatic_mutation(self):
              somatic_maf,label = calJaccardMtrx(base_path=self.somatic_mutation_maf_dir)
              label = label.map(dict(zip(
                                      sorted(label.unique()),
                                      ['peachpuff', 'lightblue']
                                      )))

              plt.figure(figsize=(7, 5))
              ax = plt.subplot(111)
              dendrogramPlot(df= 1 - somatic_maf,label=label,ax=ax)
              plt.subplots_adjust(right=.7, bottom=0.2)
              plt.savefig(
                  '{0}/somatic_mutation_cluster_tnsnv.{1}'.format(self.out_dir,self.fmt), dpi=self.DPI)
              plt.clf()

              ## Mutation profile
              with PdfPages('{}/somatic_mutation_profile_tnsnv.pdf'.format(self.out_dir)) as pdf:
                  for sampleID,tri_mtrx in iterMaf(base_path=self.somatic_mutation_maf_dir, ref_fasta=self.ref_fasta):
                      ref_mtrx = pd.concat([pd.read_table("{0}/{1}.mtrx".format(REF,x)) for x in self.cancer],axis=0)
                      tri_mtrx['Group'] = sampleID
                      fig = plot96Mtrx(df=pd.concat(
                          [tri_mtrx, ref_mtrx]), rows='Group')
                      pdf.savefig(fig,dpi=self.DPI)
                      plt.clf()
