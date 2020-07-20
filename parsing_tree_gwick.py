import pathlib
import sys
import string
import random
import os
import pandas as pd

from dollo_tree import TreeNode
import _delete as dl
import numpy as np

from graphviz import Digraph
from pandas import DataFrame
from visualize_score import plot_cell_cn_profile
import matplotlib.pyplot as plt

LIB = []


def get_random_string(length=7):
    """
    Generate a random string for naming the nodes 

    :param length: length of the string.
    :param library: all used names
    :return: A string as names for unamed nodes.
    """
    letters = string.ascii_lowercase.join("0123456789")
    result_str = ''.join(random.choice(letters) for i in range(length))
    while result_str in LIB:
        result_str = ''.join(random.choice(letters) for i in range(length))
    LIB.append(result_str)
    return result_str
    # print("Random string of length", length, "is:", result_str)


def loads(s):
    """
    Load a list of trees from a Newick formatted string.

    :param s: Newick formatted string.
    :param strip_comments: Flag signaling whether to strip comments enclosed in square \
    brackets.
    :param kw: Keyword arguments are passed through to `Node.create`.
    :return: List of Node objects.
    """
    s = s.replace(")1", ")")
    #print(s)
    return [parse_node(ss.strip()) for ss in s.split(';') if ss.strip()]


def dumps(trees):
    """
    Serialize a list of trees in Newick format.
    :param trees: List of TreeNode objects or a single Node object.
    :return: Newick formatted string.
    """
    if isinstance(trees, TreeNode):
        trees = [trees]
    return ';\n'.join([tree.newick for tree in trees]) + ';'


def load(fp):
    """
    Load a list of trees from an open Newick formatted file.

    :param fp: open file handle.
    :return: List of Node objects.
    """
    return loads(fp.read())


def read(fname, encoding='utf8'):
    """
    Load a list of trees from a Newick formatted file.
    :param fname: file path.
    :return: List of Node objects.
    """
    with pathlib.Path(fname).open(encoding=encoding) as fp:
        return load(fp)


def _parse_name_and_length(s):
    length = None
    if ':' in s:
        s, length = s.split(':', 1)
    name = get_random_string()
    return s or name, length or None


def _parse_siblings(s):
    """
    http://stackoverflow.com/a/26809037
    """
    bracket_level = 0
    current = []

    # trick to remove special-case of trailing chars
    for c in (s + ","):
        if c == "," and bracket_level == 0:
            yield parse_node("".join(current))
            current = []
        else:
            if c == "(":
                bracket_level += 1
            elif c == ")":
                bracket_level -= 1
            current.append(c)


def parse_node(s):
    """
    Parse a Newick formatted string into a `TreeNode` object.
    :param s: Newick formatted string to parse.
    :param strip_comments: Flag signaling whether to strip comments enclosed in square \
    brackets.
    :return: `TreeNode` instance.
    """
    #print(s)
    s = s.strip()
    parts = s.split(')')
    if len(parts) == 1:
        descendants, label = [], s
    else:
        if not parts[0].startswith('('):
            raise ValueError('unmatched braces %s' % parts[0][:100])
        descendants = list(_parse_siblings(')'.join(parts[:-1])[1:]))
        label = parts[-1]
    length, name = _parse_name_and_length(label)
    new_tree = TreeNode(name=name, label=length)
    new_tree.children = descendants
    print("new_tree", new_tree.label)
    return new_tree


def visualize_tree(tree, fname, cn_max, num_bins, clones = None, df=None, display=None, bins=None, chromosome=None):
    """
    visualize the tree and save the result as pdf file.
    3 modes:
        graph: visualize each node using scatter plot.
            The graph for each node will be saved in the folder "pasimony_figres/fname"
        table: visualize each node with data table
        None: visualize each node with its name
    :param tree: TreeNode object to be visualized
    :param fname: filename for the saved gv file as well as pdf
    :param cn_max: maximum copy number
    :param num_bins: total number of bins in the data
    :param df: (only used when display == 'graph') input DataFrame for visualizing the nodes into graph
    :param display: diplay mode ("graph" or "table", or "None")
    :param bins: (optional) index of bins you wish to visualize (a list of integers)
    :param chromosome: (optional for display == "graph") choose specific chromosome you want to visualize

    """
    dot = Digraph(comment='Tree for %s' % fname)
    i = 0
    if clones is None:
        vis_nodes = [node.label for node in tree.nodes]
    for node in tree.nodes:
        for child in node.children:
            if child.label in vis_nodes:
                # different modes of displaying the node
                if display == None:
                    node_display = node.label
                    child_display = child.label
                elif display == "table":
                    node_display = node_vis_table(node, num_bins, bins, cn_max)
                    child_display = node_vis_table(child, num_bins, bins, cn_max)
                elif display == "graph":
                    node_display = node_vis_graph(df, node, num_bins, bins, cn_max, chromosome, fname, s=5)
                    child_display = node_vis_graph(df, child, num_bins, bins, cn_max, chromosome, fname, s=5)
                dot.node(node.label, node_display)
                dot.node(child.label, child_display)
                dot.edge(node.label, child.label)
    dot.render('%s.gv' % fname, view=True)


def node_vis_table(node, num_bins, bins, cn_max):
    """
    gereating html code for visualizing the graph
    :param node: the node for graph generation
    :param bins: the bins we want to look at
    :param cn_max: maximum copy number
    :return: an html string
    """
    first_row = '''<tr><td>bin \ cp# </td> '''
    for i in range(cn_max + 1):
        first_row += '''<td> %d </td>''' % i
    first_row += '''</tr>'''
    html_table = '''<<table><tr><td> node: %s </td></tr>''' % node.label + first_row
    if bins == None:
        bins = [i for i in range(num_bins)]
    for bin in bins:
        one_row = '''<tr><td> %s </td>''' % bin  # scope='row'
        parsimony_scores = node.cn_score[bin, :]
        cn_backtrack = node.cn_backtrack[bin]
        # print("node", node.label)
        # print("bin_backtrack", type(cn_backtrack))
        # print("bin_backtrack1", cn_backtrack)
        for cn in range(cn_max + 1):
            score = parsimony_scores[cn]
            if cn in cn_backtrack[0].tolist():
                one_row = one_row + '''<td BGCOLOR='lightgrey'> %s </td>''' % str(score)
            else:
                one_row = one_row + '''<td> %s </td>''' % str(score)
        one_row += '''</tr>'''
        html_table += one_row
    html_table += '''</table>>'''
    # print(html_table)
    with open("file.txt", "w") as f:
        f.write(html_table)
    f.close()
    return html_table


def extracting_score(df, node, num_bins, name_tree):
    """
    extracting parsimony score as well as inferred copy number
    (as an input for node_vis_graph)
    from the dataset and save it to a csv file.
    :param df: data frame to be extracted
    :param node: input TreeNode object
    :param num_bins: total number of bins in the dataset
    :return:
    a DataFrame object containing parsimony scores and inferred copy number
    of the tree node
    """
    assert (len(node.cn_backtrack) == num_bins), "not enough bins for back trackers"
    copy_number = [node.cn_backtrack[bin][0].tolist()[0] for bin in range(num_bins)]
    # print(copy_number)
    parsimony_score = [node.cn_score[bin, cn] for (bin, cn) in enumerate(copy_number)]
    df["copy_number"] = np.array(copy_number)
    df["parsimony_score"] = np.array(parsimony_score)
    dir = file_dir(name_tree, "")
    dir = os.path.join(dir, "score_%s.csv" % node.label)

    df.to_csv(dir)
    return df


def node_vis_graph(df, node, num_bins, bins, cn_max, chromosome, name_tree, s=5):
    '''
    save and visualize current node using pyplot
    :param df: data frame for current node
    :param node: TreeNode object for current node
    :param num_bins: total number of bins
    :param cn_max: max copy number for computational efficiency
    :param chromosome: if there is a specific chromosome to visualize
    :param name_tree: name of the tree
    :param s: size of each point
    :return:
    html_image: html code for the saved image
    '''
    df = extracting_score(df, node, num_bins, name_tree).copy()
    fig = plt.figure(figsize=(16, 4))
    ax = plt.gca()
    plt.title(node.label)
    plot_cell_cn_profile(ax, df, "parsimony_score", "copy_number", cn_max, chromosome, s)
    dir = file_dir(name_tree)
    dir = os.path.join(dir, "%s.png" % node.label)
    plt.savefig(dir)
    plt.close(fig)
    html_image = '''<<TABLE>
   <TR><TD><IMG SRC=\"%s\"/></TD></TR></TABLE>>''' % dir
    return html_image


def file_dir(name_tree, sub_dir = "figure"):
    '''
    create a directory for saving figures in the current working directory
    :param name_node: name of current node
    :param name_tree: name of the tree
    :param figures: is this directory for saving figures
    :return:
    dir: directory of tree generated
    '''
    cwd = os.getcwd()

    path = os.path.join(cwd, "parsimony_%s"%sub_dir)
    if not os.path.exists(path):
        os.mkdir(path)
    path = os.path.join(path, name_tree)
    if not os.path.exists(path):
        os.mkdir(path)
    return path

def read_csv_cn(filename, nodes):
    """
    read the csv file containing copy number data
    :param filename: file name of the csv ("filename.csv")
    :param nodes: list of nodes contained in the file
    :return:
    df_location: a dataframe containing the location information of bins
    leaf_cns: a dataframe containing the copy number of each bins for all leaves
    num_bins: total number of bins in the dataset
    """
    df = pd.read_csv(filename, sep='\t', header=0)
    # data cleaning
    df['chr'] = df['chr'].astype(str)
    df = df.loc[df["start"] > 0]
    df["start"] = df["start"].astype(float)
    df["end"] = df["end"].astype(float)
    leaf_cns = {}
    for node in nodes:
        leaf_cns[node] = df[node]
    num_bins = len(df.index)
    # returning the location data only
    df_location = df.drop(columns=nodes).copy()
    print(df_location.chr.unique())
    return df_location, leaf_cns, num_bins


def sankoff_parimony(tree_fname, cp_fname,
                     segment_length,
                     cn_max, bins, clone_fname = None, display=None, vis=False):
    # read the tree from newick format
    print("---- Reading the tree : %s.newick----" % tree_fname)
    tree = read(tree_fname + ".newick")[0]
    nodes = [node.label for node in tree.leaves]

    # read the copy number data in specified file
    print("---- Reading copy number data: %s ----" % cp_fname)
    df_location, leaf_cns, num_bins = read_csv_cn(cp_fname, nodes)

    # calculate the parsimony score
    print("---- Calculating parsimony score: %s ----" % tree_fname)
    dl.calc_score_recursive_vect(tree, leaf_cns,
                                 segment_length,
                                 cn_max, num_bins)
    dl.cn_backtrack_root(tree, num_bins)
    if clone_fname is not None:
        print("---- Labeling Clones using: %s ----" % clone_fname)
        _, clones = label_clones(clone_fname, tree)


    if vis:
        # visualize the tree
        print("---- Visualizing the tree: %s ----" % tree_fname)
        # pt.visualize_tree(tree, tree_fname, cn_max, num_bins, display, bins)
        if clone_fname is not None:
            visualize_tree(tree, tree_fname, cn_max, num_bins, clones, df_location, display, bins)
        else:
            visualize_tree(tree, tree_fname, cn_max, num_bins, df_location, display, bins)

    # save the score of the root in a csv file:
    print("---- Saving the score: %s ----" % tree_fname)
    pd.DataFrame(tree.cn_score.transpose(), columns=[str(i + 1) for i in range(num_bins)]).to_csv(tree_fname + ".csv")

def read_assignment(fname):
    """
    read the clone assignment files
    :param fname: filename for clone assignments
    :return:
    assignments: the dictionary object with key as clone names,

    clones: list of clones

    """
    df = pd.read_csv(fname, sep='\t', header=0)
    clones = df['clone_id'].unique().tolist()
    clones.remove("None")
    assignment = {}
    for clone in clones:
        cells = df.loc[df["clone_id"] == clone]["cell_id"]
        assignment[clone] = set(cells.tolist())
    return assignment, clones

def label_clones(fname, tree):
    assignment, clones = read_assignment(fname)
    num_cluster = len(assignment)
    #print("num_clusters", num_cluster)
    #print("cluster", clones)
    bfs([], [], tree, num_cluster, assignment)
    return assignment, clones

def bfs(queue, visited, tree, count, assignment):
    clusters = assignment.values()
    #print(clusters)
    visited.append(tree)
    queue.append(tree)

    while (count != 0):
        s = queue.pop(0)
        #print("current node", s.label, count)
        #print("que", [node.label for node in queue])
        #print("visited", [node.label for node in visited])
        queue_name = [node.label for node in queue]
        visited_name = [node.label for node in visited]
        for child in s.children:
            if child.label not in visited_name:
                visited.append(child)
                queue.append(child)
        leave_names = set([leaf.label for leaf in s.leaves])
        #print("leaves", leave_names)
        if leave_names in clusters:
            print("Yes")
            s.label = get_key(leave_names, assignment)
            count -= 1
    #print("total visits", len(visited))

def get_key(val, my_dict):
    for key, value in my_dict.items():
         if val == value:
             return key

def main():
    # num_site = int(sys.argv[1])
    # cn_max = int(sys.argv[2])
    '''
    tree_gw = "((A,G), ((C,(D,B)1:0),(F,E)1:0)1:0.01)"
    tree1 = loads(tree_gw)[0]
    assignment = {"X": set(["A","G"]), "Y":set(["B", "C", "D"]), "Z": set(["F","E"])}
    tree1.print_attr_map("name")
    #label_clones(assignment, tree1)
    #tree1.print_attr_map("name")
    '''
    '''
    read_assignment("SA1053_clones.tsv")
    tree1 = read("SA1053_sitka_hdbscan_tree.newick")[0]
    label_clones("SA1053_clones.tsv", tree1)
    print(tree1.label)
    clones = ["A", "B", "C", "D", "E", "F", "G", "H"]
    print(tree1.newick_str("name"))
    for node in tree1.nodes:
        if node.label in clones:
            print(node.label)
            node.count_leaves()
            print(node.num_leaves)
    '''

    cp_fname = "SA1053_filtered_cell_states.tsv"
    tree_name = "SA1053_sitka_hdbscan_tree"#"SA1053_clone_tree"
    clone_fname = "SA1053_clones.tsv"
    segment_length = 1
    cn_max = 12
    bins = [0, 500, 600, 700, 800]
    display = "graph"

    # print(signature(pt.visualize_tree))
    sankoff_parimony(tree_name, cp_fname,
                     segment_length,
                     cn_max, bins,
                     clone_fname,
                     display, vis=True)

    '''
    tree_gw = "((E,D),((C,(A,B)),(F,G)))"
    tree1 = loads(tree_gw)[0]
    print("Tree1")
    cn_max = 4
    num_site = 2
    leaf_cns = {}
    tree1.count_leaves()

    for node in tree1.leaves:
        leaf_cns[node.label] = np.random.randint(cn_max, size=num_site)

    dl.calc_score_recursive_vect(tree1, leaf_cns, 1, cn_max, num_site)
    print(tree1.cn_score)
    print("number of leaf: %d" % tree1.num_leaves)
    print("number of sites: %d" % num_site)
    print("number of max copies: %d" % cn_max)
    tree1.cn_backtrack = node.cn_backtrack

    visualize_tree(tree1, "name", cn_max, [0, 1], "graph")


    for node in tree1.nodes:
        print("name", node.label)
        print("backtrack", node.cn_backtrack)
    '''


main()
