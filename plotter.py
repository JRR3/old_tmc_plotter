#########################################################
#Princess Margaret Cancer Research Tower
#Schwartz Lab
#Javier Ruiz Ramirez
#July 2024
#########################################################
#This is a Python script to produce TMC trees using
#the original too-many-cells tool.
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7439807/
#########################################################
#Questions? Email me at: javier.ruizramirez@uhn.ca
#########################################################
import subprocess
import pandas as pd
import matplotlib as mpl
import os

class TMC:
    def __init__(self):
        #High-High, High-Low, Low-High, Low-Low.
        self.high_low_colors = ['purple',
                                'red',
                                'blue',
                                'aqua']
        self.high_low_colors_hex = []
        self.list_of_genes = ["FAP","TCF21"]

        self.list_of_colors = ["red", "blue"]

        self.gene_to_tuple = {}

        self.matrix_path = "./tmci_mtx_data"

        self.tmc_tree_path = "./data_to_regenerate_tree"

        self.method = "MadMedian"

        self.threshold = 1.5

        self.all_genes = []

        self.draw_node_numbers = False

        self.cell_annotations_fname = "cell_labels.csv"

        self.output = "."
        self.output = os.path.join(self.output, "tmc_out")

        self.saturation = 2

        os.makedirs(self.output, exist_ok=True)

        self.tree_output = os.path.join(self.output, "tree.svg")
        self.tree_output = "tree.svg"

        self.cluster_path = os.path.join(self.output, "clusters.csv")


    def create_gene_objects(self):

        for gene, color in zip(self.list_of_genes, self.list_of_colors):
            color_hex = mpl.colors.cnames[color]
            color = '\\\"' + color_hex.lower() + '\\\"'
            gene_mod = '\\\"' + gene + '\\\"'
            self.gene_to_tuple[gene] = (gene_mod, color)


        for gene, T in self.gene_to_tuple.items():
            gene_label, hex_color = T
            gene_bicolor = '(' + gene_label + ', ' 
            gene_bicolor += self.method
            gene_bicolor += ' '
            gene_bicolor += str(self.threshold)
            gene_bicolor += ')'
            self.all_genes.append(gene_bicolor)

        for color in self.high_low_colors:
            color_hex = mpl.colors.cnames[color]
            color_hex = '\\\"' + color_hex.lower() + '\\\"'
            self.high_low_colors_hex.append(color_hex)

        self.color_list = '[' + ','.join(self.high_low_colors_hex) + ']'
        self.gene_list = '[' + ','.join(self.all_genes) + ']'
        #'\"DrawItem (DrawThresholdContinuous [(\\\"FH\\\", Exact 0), (\\\"FL\\\", Exact 0)])\"',
        self.gene_txt = '\"DrawItem (DrawThresholdContinuous '
        self.gene_txt += self.gene_list
        self.gene_txt += ')\"'

    def execute_command(self):

        print(self.gene_txt)

        if self.draw_node_numbers:
            node_numbers = '--draw-node-number'
        else:
            node_numbers = ''
        command = ['too-many-cells',
                'make-tree',
                '--matrix-path',
                self.matrix_path,
                '--prior',
                self.tmc_tree_path,
                #'--normalization',
                #'UQNorm'
                '--labels-file',
                self.cell_annotations_fname,
                node_numbers,
                '--draw-mark',
                'MarkModularity',
                '--feature-column',
                '1',
                '--draw-leaf',
                #'\"DrawItem (DrawThresholdContinuous [(\\\"FH\\\", Exact 0), (\\\"FL\\\", Exact 0)])\"',
                self.gene_txt,
                '--draw-colors',
                self.color_list,
                '--draw-scale-saturation',
                str(self.saturation),
                '--dendrogram-output',
                self.tree_output,
                '--labels-output',
                '--output',
                self.output,
                '>',
                self.cluster_path]
        command = list(filter(len, command))
        command = ' '.join(command)
        #print(">>>>>>>>>>>")
        #print(command)
        p = subprocess.call(command, shell=True)

    def run_tmc(self):
        self.create_gene_objects()
        self.execute_command()

obj = TMC()
obj.run_tmc()
