import argparse
import sys
from random import randint
from math import sqrt, log
from copy import copy
import networkx as nx
import numpy as np
import operator

reads = {} # reads = {read1: [protein1, protein2, ...]} hashmap storing the list of proteins to which a read maps
proteins = {} # proteins = {protein1: [read1, read2, ...]} hashmap storing the list of reads that mapped to it
path_expression = {} # path_expression = {enzyme: exp} "expression" of enzymes
new_path_expression = {} # this will store the same thing as path_expression, but updated after each iteration
names_map = {}
weights_map = {}
expression = {} #stores enzyme ec_expression
new_expression = {} #stores enzymes expression after each iteration
new_expression_out = 0

INF = 10000000

def parse_arguments():
    parser = argparse.ArgumentParser(description='EMPathways2 - A tool for pathway reconstruction from metagenomic data.')
    parser.add_argument('-ge', '--gene_enzyme_file', required=True, help='Gene-enzyme relation file') #gene-enzyme-weight file with 1's
    parser.add_argument('-epd', '--enzyme_pathway_dictionary', required=True, help='Enzyme Pathway Dictionary (EPD) file') #enzyme-pathway dictionary file
    parser.add_argument('-ge', '--gene_expression_file', required=True, help='Gene expression file') #gene_fpkm IsoEM output file
    parser.add_argument('-eo', '--enzyme_output', required=False, help='Filename to write Enzymes to')
    parser.add_argument('-po', '--pathway_output', required=False, help='Filename to write Pathways to')
    parser.add_argument('--theta_e', type=float, default=0.000001, help='Theta for enzyme convergence')
    parser.add_argument('--theta_p', type=float, default=0.000001, help='Theta for pathway convergence')
    return parser.parse_args()

args = parse_arguments()
ec_pathway_dict = args.enzyme_pathway_dictionary
ge_file = args.gene_enzyme_file
weights_file = args.weights_file
outputenzymes = args.enzyme_output
outputpathways = args.pathway_output
theta_e = args.theta_e
theta_p = args.theta_p
numSample = weights_file[-7:-5]

print("Sample:", numSample)

#Utility Functions
def error(new_expression, expression):
    error = 0
    for key in expression:
        if (error < abs(expression[key] - new_expression[key])):
            error = abs(expression[key] - new_expression[key])
    return error

def error_r(r_ij_new, r_ij):
    error = 0
    for ec, pathways in r_ij.items():
        for pathway in pathways:
            if (error < abs(r_ij_new[ec][pathway] - r_ij[ec][pathway])):
                error = abs(r_ij[ec][pathway] - r_ij_new[ec][pathway])
    return error

def load_ec_pathway_dict(ec_pathway_dict):
    with open(ec_pathway_dict) as f:
        pd = f.readlines()

    pd = list(map(lambda x: x.strip().split("?"), pd))

    pathways = []
    all_enzymes = []

    pathway_ec_dict = {}
    ec_pathway_dict = {}
    enzyme_list = {}

    for x, y in pd:
        paths = list(set(y.split(",")))
        enzyme = x
        ec_pathway_dict[x] = paths
        for path in paths:
            if path not in pathway_ec_dict:
                pathway_ec_dict[path] = []
            pathway_ec_dict[path].append(x)
            pathways.append(path)
            if enzyme not in all_enzymes:
                all_enzymes.append(enzyme)

    return ec_pathway_dict, pathway_ec_dict, pathways, all_enzymes

def load_weights_file(weights_file):
    with open(weights_file) as f:
        weights = f.readlines()
        weights = list(map(lambda x: x.strip().split(), weights))

    weights_map = {}
    for x, y in weights:
        weights_map[x] = float(y)

    return weights_map

def load_alignment_file(ge_file, ec_pathway_dict, read_dict, enzyme_dict):
    edges = []
    weights = []
    ecnumbers_global = []

    with open(ge_file) as f:
        for line in f.readlines():
            line = line.split()
            read, ecnumber = line[0], line[1]
            if ecnumber not in ec_pathway_dict:
                continue
            bitscore = float(line[2])

            if read not in read_dict:
                read_dict[read] = len(read_dict)
            if ecnumber not in enzyme_dict:
                enzyme_dict[ecnumber] = len(enzyme_dict)

            edges.append((read_dict[read], enzyme_dict[ecnumber]))
            weights.append(bitscore)
            ecnumbers_global.append(ecnumber)

    return edges, weights, ecnumbers_global


def construct_weight_matrix(edges, weights, read_dict, enzyme_dict):
    weight_matrix = np.zeros((len(read_dict), len(enzyme_dict)))

    for (read, enzyme), weight in zip(edges, weights):
        weight_matrix[read, enzyme] = weight

    return weight_matrix

def normalize_weights(weight_matrix):
    row_sums = weight_matrix.sum(axis=1)
    normalized_weight_matrix = weight_matrix / row_sums[:, np.newaxis]
    return normalized_weight_matrix.reshape((-1, 1))

def em_enzymes(weight_matrix, enzyme_count, theta_e):
    enzyme_expression = np.full(enzyme_count, 1 / enzyme_count)

    while True:
        print("Enzyme Expression:", enzyme_expression)
        new_enzyme_expression = em1_step(weight_matrix, enzyme_expression, weight_matrix.T)
        new_enzyme_expression /= np.sum(new_enzyme_expression)

        if np.max(np.abs(enzyme_expression - new_enzyme_expression)) < theta_e:
            break

        enzyme_expression = new_enzyme_expression

    return enzyme_expression


def em_pathways(weight_matrix, enzyme_pathway_matrix, pathway_count, theta_p):
    pathway_expression = np.full(pathway_count, 1 / pathway_count)

    while True:
        new_pathway_expression = em1_step(enzyme_pathway_matrix, pathway_expression, weight_matrix)
        new_pathway_expression /= np.sum(new_pathway_expression)
        if np.linalg.norm(pathway_expression - new_pathway_expression) < theta_p:
            break
        pathway_expression = new_pathway_expression

    return pathway_expression

def em_pathway_enzymes(weight_matrix, enzyme_pathway_matrix, pathway_count, enzyme_count, theta_e, theta_p, enzyme_dict, pathway_ec_dict):
    enzyme_expression = np.full(enzyme_count, 1 / enzyme_count)
    pathway_expression = np.full(pathway_count, 1 / pathway_count)

    while True:
        new_enzyme_expression, new_pathway_expression = em_two_steps(weight_matrix, enzyme_expression, enzyme_pathway_matrix, theta_e)

        if np.linalg.norm(enzyme_expression - new_enzyme_expression) < theta_e and np.linalg.norm(pathway_expression - new_pathway_expression) < theta_p:
            break

        enzyme_expression = new_enzyme_expression
        pathway_expression = new_pathway_expression

    return enzyme_expression, pathway_expression

#H = r_ij_new[enzyme][pathway]
#G = 1/len(pathways)
#F = f_i_exp[enzyme]

def em1_step(H, G, F):
    E = ((H*G).T*F).T
    E = (E.T/E.sum(axis=1)).T
    return E.sum(axis=0)/E.sum()

def em2_step(H, G, F):
    F_exp = np.matmul(H, G)
    F_exp = F_exp/F_exp.sum()
    return (H.T*F/F_exp).T

def em_two_steps(H, G, F, e):
    #EM1
    while True:
        G_new = em1_step(H, G, F)
        if np.linalg.norm(G-G_new) < e:
            break
        G = G_new
    #EM2
    while True:
        H_new = em2_step(H, G, F)
        if np.linalg.norm(H-H_new) < e:
            break
        H = H_new
    return G_new, H_new


def main():
    
    read_dict = {}
    enzyme_dict = {}
    iteration = 0

    ec_pathway_dict, pathway_ec_dict, pathways, all_enzymes = load_ec_pathway_dict(args.enzyme_pathway_dictionary)
    edges, weights, ecnumbers_global = load_alignment_file(args.gene_enzyme_file, ec_pathway_dict, read_dict, enzyme_dict)
    weight_matrix = construct_weight_matrix(edges, weights, read_dict, enzyme_dict)
    normalized_weight_matrix = normalize_weights(weight_matrix)
    
    print("Iteration:", iteration)
    enzyme_expression = em_enzymes(normalized_weight_matrix, len(enzyme_dict), theta_e)

    enzyme_pathway_matrix = np.zeros((len(enzyme_dict), len(pathways)))

    enzyme_pathway_matrix = np.zeros((len(enzyme_dict), len(pathways)))
    for enzyme, pathways in ec_pathway_dict.items():
        if enzyme in enzyme_dict:
            enzyme_idx = enzyme_dict[enzyme]
            for pathway in pathways:
                pathway_idx = pathways.index(pathway)
                enzyme_pathway_matrix[enzyme_idx, pathway_idx] = 1

    enzyme_pathway_matrix = enzyme_pathway_matrix.T
    pathway_expression = em_pathways(normalized_weight_matrix, enzyme_pathway_matrix, len(pathway_ec_dict), theta_p)

    
    # Run em_pathway_enzymes() until convergence
    while True:
        new_enzyme_expression, new_pathway_expression = em_pathway_enzymes(normalized_weight_matrix, enzyme_pathway_matrix, len(pathway_ec_dict), len(enzyme_dict), theta_e, theta_p, enzyme_dict, pathway_ec_dict)

        if error(new_enzyme_expression, enzyme_expression) < theta_e and error(new_pathway_expression, pathway_expression) < theta_p:
            break

        enzyme_expression = new_enzyme_expression
        pathway_expression = new_pathway_expression

        print("Enzyme expression:")
        for enzyme, expression_value in zip(enzyme_dict.keys(), enzyme_expression):
            print(enzyme, expression_value)

        print("Pathway expression:")
        for pathway, expression_value in zip(pathway_ec_dict.keys(), pathway_expression):
            print(pathway, expression_value)

    # Write enzyme expressions to output file, if provided
    if outputenzymes:
        with open(outputenzymes, 'w') as f:
            f.write("Enzyme expression:\n")
            for enzyme, expression_value in zip(enzyme_dict.keys(), enzyme_expression):
                f.write(enzyme + " " + str(expression_value) + "\n")

    # Write pathway expressions to output file, if provided
    if outputpathways:
        with open(outputpathways, 'w') as f:
            f.write("Pathway expression:\n")
            for pathway, expression_value in zip(pathway_ec_dict.keys(), pathway_expression):
                f.write(pathway + " " + str(expression_value) + "\n")

if __name__ == '__main__':
    main()

       
