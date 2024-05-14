from utils.viziumHD import viziumHD

def run(path,pathout):
    viziumHD(path = path,outPath = pathout,totalThr = 1500,geneThr = 1500 , bins_gene_plot = 20,bins_gene_plot_thr = 20,qcFilePrefix = 'qcFilePrefix')
    
if __name__ == "__main__":
    path = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Mouse_Brain_square_example/square_008um"
    pathout = "/data/kanferg/Sptial_Omics/playGround/out"
    run(path,pathout)
    

'''
self.totalThr = totalThr
self.bins_total = bins_total
self.bins_gene_plot = 200
self.geneThr = geneThr
self.bins_gene_plot_thr = bins_gene_plot_thr
self.qcFilePrefix = qcFilePrefix
'''