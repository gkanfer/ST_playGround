import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import os
import gzip
import numpy as np

class viziumHD:
    '''
    A class to handle the processing and quality control reporting of Visium high-definition spatial transcriptomics data.

    Parameters
    ----------
    path : str
        The path to the directory containing the spatial transcriptomics data files.
    outPath : str
        The output directory path where the results, including the quality control report, will be saved.
    totalThr : int, optional
        Threshold for filtering out cells based on total counts below this threshold. Default is 10000.
    bins_total : int, optional
        Number of bins to use for the histogram of total counts. Default is 40.
    bins_gene_plot : int, optional
        Number of bins to use for the histogram plotting the number of genes by counts. Default is 60.
    geneThr : int, optional
        Threshold for filtering out cells based on the number of genes detected below this threshold. Default is 4000.
    bins_gene_plot_thr : int, optional
        Number of bins to use for the histogram of the number of genes by counts with a threshold. Default is 60.
    '''
    def __init__(self,path,outPath,totalThr = 10000, bins_total = 40, bins_gene_plot = 60, geneThr = 4000,
                 bins_gene_plot_thr = 60, qcFilePrefix = ""):
        self.path = path
        self.outPath = outPath
        self.totalThr = totalThr
        self.bins_total = bins_total
        self.bins_gene_plot = bins_gene_plot
        self.geneThr = geneThr
        self.bins_gene_plot_thr = bins_gene_plot_thr
        self.qcFilePrefix = qcFilePrefix
        self.filterNorm = filterNorm
        self.parquet_to_csv()
        self.andata = self.readVizHD()
        self.qcReport()
        
    def parquet_to_csv(self):
        file_path = os.path.join(self.path,'spatial/tissue_positions_list.csv')
        # Read the Parquet file
        if os.path.exists(file_path):
            return
        else:
            df = pd.read_parquet(os.path.join(self.path,'spatial/tissue_positions.parquet'))
            # Write to a CSV file
            df.to_csv(os.path.join(self.path,'spatial/tissue_positions_list.csv'), index=False)
    
    def readVizHD(self):
        return sc.read_visium(path = self.path)
    
    def filterANDnorm(self,**kwargs):
        sc.pp.filter_cells(self.andata,**kwargs)
        sc.pp.filter_genes(self.andata,**kwargs)
        print("normalize total")
        sc.pp.normalize_total(self.andata)
        print("log transform")
        sc.pp.log1p(self.andata)
        print("scale")
        sc.pp.scale(self.andata,**kwargs)  
    
    def qcReport(self):
        sc.pp.calculate_qc_metrics(self.andata, inplace=True)
        print("start qc report")
        with PdfPages(os.path.join(self.outPath, f'Quality_Control_{self.qcFilePrefix}.pdf')) as pdf:
            plt.rcParams['figure.dpi'] = 150
            plt.rcParams['font.family'] = ['serif']
            plt.rcParams['font.size'] = 12
            plt.rcParams['axes.labelsize'] = 12
            plt.rcParams['axes.titlesize'] = 12
            plt.rcParams['xtick.labelsize'] = 12
            plt.rcParams['ytick.labelsize'] = 12
            fig, axs = plt.subplots(1, 4, figsize=(20, 5))  # Adjusted figsize for better readability

            # Plot for total counts
            sns.histplot(self.andata.obs["total_counts"], kde=False, ax=axs[0])
            axs[0].set_title('Total Counts per Cell')
            axs[0].set_xlabel('Total Counts')
            axs[0].set_ylabel('Frequency')

            # Plot for total counts with a threshold
            sns.histplot(
                self.andata.obs["total_counts"][self.andata.obs["total_counts"] < self.totalThr],
                kde=False,
                bins=self.bins_total,
                ax=axs[1],
            )
            axs[1].set_title(f'Total Counts per Cell (Threshold < {str(self.totalThr)})')
            axs[1].set_xlabel('Total Counts')
            axs[1].set_ylabel('Frequency')

            # Plot for number of genes by counts
            sns.histplot(self.andata.obs["n_genes_by_counts"], kde=False, bins=self.bins_gene_plot, ax=axs[2])
            axs[2].set_title('Number of Genes Detected per Cell')
            axs[2].set_xlabel('Number of Genes')
            axs[2].set_ylabel('Frequency')

            # Plot for number of genes by counts with a threshold
            sns.histplot(
                self.andata.obs["n_genes_by_counts"][self.andata.obs["n_genes_by_counts"] < self.geneThr],
                kde=False,
                bins=self.bins_gene_plot_thr,
                ax=axs[3],
            )
            axs[3].set_title(f'Number of Genes Detected per Cell (Threshold < {str(self.geneThr)})')
            axs[3].set_xlabel('Number of Genes')
            axs[3].set_ylabel('Frequency')

            fig.tight_layout()
            pdf.savefig()
            plt.close()
        
        
             