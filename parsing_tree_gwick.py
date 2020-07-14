import pathlib
import sys
import string
import random

from dollo_tree import TreeNode
import _delete as dl
import numpy as np

from graphviz import Digraph

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
    print(s)
    s = s.strip()
    parts = s.split(')')
    if len(parts) == 1:
        descendants, label = [], s
    else:
        if not parts[0].startswith('('):
            raise ValueError('unmatched braces %s' % parts[0][:100])
        descendants = list(_parse_siblings(')'.join(parts[:-1])[1:]))
        label = parts[-1]
    name, length = _parse_name_and_length(label)
    new_tree = TreeNode(name=name, label=length)
    new_tree.children = descendants
    # print("new_tree", new_tree.name)
    return new_tree


def visualize_tree(tree, fname, cn_max, bin=None):
    """

    :param tree:
    :param fname:
    :param bin:
    :return:
    """
    dot = Digraph(comment='Tree for %s' % fname)
    i = 0
    for node in tree.nodes:
        for child in node.children:
            node_html = node_vis(node, bin, cn_max)
            child_html = node_vis(child, bin, cn_max)
            dot.node(node.name, node_html)
            dot.node(child.name, child_html)
            dot.edge(node.name, child.name)
    dot.view()


def node_vis(node, bins, cn_max):
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
    html_table = '''<<table><tr><td> node: %s </td></tr>''' % node.name + first_row
    print("bins", bins)
    for bin in bins:
        one_row = '''<tr><td> %s </td>''' % bin  # scope='row'
        parsimony_scores = node.cn_score[bin, :]
        cn_backtrack = node.cn_backtrack[bin]
        print("node", node.name)
        print("bin_backtrack", type(cn_backtrack))
        print("bin_backtrack1", cn_backtrack)
        for cn in range(cn_max + 1):
            score = parsimony_scores[cn]
            if cn in cn_backtrack[0].tolist():
                one_row = one_row + '''<td BGCOLOR='lightgrey'> %s </td>''' % str(score)
            else:
                one_row = one_row + '''<td> %s </td>''' % str(score)
        one_row += '''</tr>'''
        html_table += one_row
    html_table += '''</table>>'''
    print(html_table)
    with open("file.txt", "w") as f:
        f.write(html_table)
    f.close()
    return html_table


def main():
    # num_site = int(sys.argv[1])
    # cn_max = int(sys.argv[2])
    #fname = "OV2295-cn-tree.newick"
    fname = "/work/shah/funnellt/projects/sc-mutsig/analysis/corrupt_tree/2295/sitka_hdbscan_tree.newick"
    #tree_gw = "(a,(b, c, d))"
    tree1 = read(fname)[0]
    tree1 = loads(tree_gw)[0]
    print("Tree1")
    # tree1.print_attr_map("name")
    # tree2 = loads(tree)[0]
    # print(0 + np.inf)
    cn_max = 4
    num_site = 2
    leaf_cns = {}
    tree1.count_leaves()
    for node in tree1.leaves:
        leaf_cns[node.name] = np.random.randint(cn_max, size=num_site)

    # print(leaf_cns)
    # leaf_cns = {"a":[1,2], "b":[2,4], "c":[3,6], "d":[4,8]}
    # dl.calc_score_recursive(tree1, leaf_cns1, 1, 8)
    # print(tree1.cn_score)
    dl.calc_score_recursive_vect(tree1, leaf_cns, 1, cn_max, num_site)
    print(tree1.cn_score)
    print("number of leaf: %d" % tree1.num_leaves)
    print("number of sites: %d" % num_site)
    print("number of max copies: %d" % cn_max)
    tree1.cn_backtrack = node.cn_backtrack
    '''
    for node in tree1.nodes:
        print("name", node.name)
        print("backtrack", node.cn_backtrack)
    '''
    visualize_tree(tree1, "name", cn_max, [0, 1])


main()
