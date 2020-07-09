import itertools
import json
import pathlib

class TreeNode(object):
    def __init__(self, name=None, label=0):
        self.children = []
        self.name = name
        self.label = label
        self.num_leaves = 0


    def __eq__(self, other):
        return self.label == other.label

    def __contains__(self, other):
        for node in self.nodes:
            if node == other:
                return True

        return False

    def copy(self):
        node = TreeNode(self.name, self.label)
        for child in self.children:
            node.children.append(child.copy())
        return node

    @property
    def is_leaf(self):
        return len(self.children) == 0

    @property
    def leaves(self):
        '''
        List of leaves in the tree rooted at self.
        '''
        if len(self.children) == 0:
            yield self

        else:
            for child in self.children:
                for leaf in child.leaves:
                    yield leaf

    def count_leaves(self):
        '''
        List of leaves in the tree rooted at self.
        '''
        if len(self.children) == 0:
            self.num_leaves = 1

        else:
            for child in self.children:
                for leaf in child.leaves:
                    leaf.count_leaves()
                    self.num_leaves = self.num_leaves + leaf.num_leaves
                    print(leaf.name, self.num_leaves)

    @property
    def nodes(self):
        '''
        List of nodes in the tree rooted at self.
        '''
        yield self
        for descendent in self.descendents:
            yield descendent


    @property
    def descendents(self):
        '''
        List of descendent nodes in the tree rooted at self.
        '''
        for child in self.children:
            for node in child.nodes:
                yield node

    @property
    def edges(self):
        for node in self.nodes:
            for child in node.children:
                yield node, child

    def relabel(self):
        '''
        Relabel each node in the tree with an integer label.
        '''
        self._relabel(itertools.count())

    def _relabel(self, labeller):
        self.label = next(labeller)

        for child in self.children:
            child._relabel(labeller)

    def find_parent(self, other_child):
        '''
        Find parent of node by label.
        '''
        for node in self.nodes:
            for child in node.children:
                if child.label == other_child.label:
                    return node

    def add_parent(self, parent=None):
        '''
        Add parent node to each node.
        '''
        self.parent = parent

        for child in self.children:
            child.add_parent(self)

    def add_ancestors(self, ancestors=()):
        '''
        Add ancestor list to each node.
        '''
        self.ancestors = ancestors

        for child in self.children:
            child.add_ancestors(ancestors + (self,))

    @property
    def leaf_name_str(self):
        if len(self.children) == 0:
            return str(self.name)

        else:
            return '(' + ','.join([child.leaf_name_str for child in self.children]) + ')'

    @property
    def label_str(self):
        if len(self.children) == 0:
            return str(self.label) + ':' + str(self.name)

        else:
            return str(self.label) + ':(' + ','.join([child.label_str for child in self.children]) + ')'

    def _newick_str(self, brlen_attr_name=None):
        if self.is_leaf:
            newick_str = ' '.join([str(self.label), str(self.name)])

        else:
            newick_str = '(' + ','.join([child._newick_str(brlen_attr_name)
                                         for child in self.children]) + ')' + str(self.label)

        if brlen_attr_name is not None:
            newick_str += ':{0}'.format(getattr(self, brlen_attr_name))

        return newick_str

    def newick_str(self, brlen_attr_name=None):
        return '(' + self._newick_str(brlen_attr_name) + ');'

    def get_attr_map(self, attr):
        attr_map = {self.label: getattr(self, attr, None)}

        for child in self.children:
            attr_map.update(child.get_attr_map(attr))

        return attr_map

    def print_attr_map(self, attr, prefix=''):
        print('{0}{1}:{2}'.format(prefix, self.label, getattr(self, attr, None)))

        for child in self.children:
            child.print_attr_map(attr, prefix + '  ')


def create_subtree(node_string):
    node = TreeNode()

    cnt = 0

    start_idx = 0

    for idx in range(len(node_string)):
        if node_string[idx] == '(':
            cnt += 1

        else:
            cnt -= 1

        if cnt == 0:
            node.children.append(create_subtree(node_string[start_idx + 1:idx]))

            start_idx = idx + 1

    return node


def augment_with_leaf(tree, leaf_name):
    for node in tree.nodes:
        internal = TreeNode()

        internal.children.append(node.copy())

        internal.children.append(TreeNode(leaf_name))

        augmented_tree = tree.copy()

        parent_node = augmented_tree.find_parent(node)

        if parent_node is None:
            augmented_tree = internal

        else:
            new_children = list()

            new_children.append(internal)

            for child in parent_node.children:
                if child.label != node.label:
                    new_children.append(child)

            parent_node.children = new_children

        augmented_tree.relabel()

        yield augmented_tree

def load(fp, strip_comments=False, **kw):
    """
    Load a list of trees from an open Newick formatted file.

    :param fp: open file handle.
    :param strip_comments: Flag signaling whether to strip comments enclosed in square \
    brackets.
    :param kw: Keyword arguments are passed through to `Node.create`.
    :return: List of Node objects.
    """
    kw['strip_comments'] = strip_comments
    return loads(fp.read(), **kw)


def has_leaf_name_groups(tree, leaf_name_groups):
    leaf_names = frozenset([leaf.name for leaf in tree.leaves])

    node_groups = set([frozenset([])])
    for node in tree.nodes:
        node_groups.add(frozenset([leaf.name for leaf in node.leaves]))

    for group in leaf_name_groups:
        if frozenset(group).intersection(leaf_names) not in node_groups:
            return False

    return True


def enumerate_labeled_trees(leaf_names, leaf_name_groups=None):
    trees = [TreeNode(leaf_names[0])]

    for leaf_idx in range(1, len(leaf_names)):
        augmented = list()

        for prev_tree in trees:
            for tree in augment_with_leaf(prev_tree, leaf_names[leaf_idx]):
                if leaf_name_groups is None or has_leaf_name_groups(tree, leaf_name_groups):
                    augmented.append(tree)

        trees = augmented

    if len(trees) == 0:
        raise ValueError('no trees generated')

    return trees



