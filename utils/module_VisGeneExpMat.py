import brewer2mpl
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


def map_color(l_stage, palette_used="Set2"):
    l_stage_uniq = []
    for stage in l_stage:
        if stage not in l_stage_uniq:
            l_stage_uniq.append(stage)

    l_stage_uniqIdx = range(len(l_stage_uniq))
    l_color = sns.color_palette(palette_used, len(l_stage_uniqIdx))
    M_color = dict(zip(l_stage_uniq, l_color))
    l_color_used = [ M_color[stage] for stage in l_stage ]
    return l_color_used

class PlotGivenGene(object):
    def __init__(self, infile_exp, infile_sam):
        self.infile_exp = infile_exp
        self.infile_sam = infile_sam
        self.__load_info()

    def __load_info(self):
        self.exp_info_mat = pd.read_csv(self.infile_exp, sep="\t", index_col=[0])
        self.sam_info = pd.read_csv(self.infile_sam, sep="\t")
    
    
    def __select_gene_from_pd(self, gene, scale_log):
        pd_exp_plot = pd.DataFrame({
            "Sample": self.sam_info['rename'].values,
            "Stage" : self.sam_info['stage'].values,
            "exp"   : self.exp_info_mat.loc[gene].values    
        })
        if scale_log:
            pd_exp_plot["exp"] = np.log10(pd_exp_plot["exp"] + 1)

        return pd_exp_plot
    
    def select_gene(self, gene, scale_log=False):
        pd_exp_plot = self.__select_gene_from_pd(gene, scale_log)
        l_color_used = map_color(pd_exp_plot['Stage'].values)
        plot_group1 = sns.barplot(x = pd_exp_plot['Sample'], y = pd_exp_plot['exp'], colorgiven=l_color_used)
        for item in plot_group1.get_xticklabels():
            item.set_rotation(45)
            item.set_ha("right")
            
        return pd_exp_plot, plot_group1

    def select_gene_groupBox(self, gene, scale_log=False):
        pd_exp_plot = self.__select_gene_from_pd(gene, scale_log)
        plot_group1 = sns.boxplot(x = pd_exp_plot['Stage'], y = pd_exp_plot['exp'])            
        return pd_exp_plot, plot_group1


    def give_pd(self, pd_given, scale_log=False):
        pd_exp_plot = pd_given
        l_color_used = map_color(pd_exp_plot['Stage'].values)
        if scale_log:
            pd_exp_plot["exp"] = np.log10(pd_exp_plot["exp"] + 1)
        
        plot_group1 = sns.barplot(x = pd_exp_plot['Sample'], y = pd_exp_plot['exp'], colorgiven=l_color_used)
        for item in plot_group1.get_xticklabels():
            item.set_rotation(45)
            item.set_ha("right")
        
        return pd_exp_plot, plot_group1


    def give_pd_groupBox(self, pd_given, scale_log=False):
        pd_exp_plot = pd_given
        plot_group1 = sns.boxplot(x = pd_exp_plot['Stage'], y = pd_exp_plot['exp'])            
        return pd_exp_plot, plot_group1
