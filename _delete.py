import gzip
import numpy as np
import pandas as pd
import pickle



import dollo_tree


def calc_score_recursive(node, leaf_cns, segment_length, cn_max):
    if node.is_leaf:
        node.cn_score = np.repeat(np.inf, cn_max + 1)

        node.cn_score[leaf_cns[node.name]] = 0.

    else:
        node.cn_score = np.zeros(cn_max + 1)

        for child in node.children:
            calc_score_recursive(child, leaf_cns, segment_length, cn_max)

            child.cn_backtrack = [[]] * (cn_max + 1)

        for cn in range(0, cn_max + 1):
            for child in node.children:
                child_cn_scores = np.repeat(0., cn_max + 1)

                for child_cn in range(0, cn_max + 1):
                    transition_score = calculate_transition_score(cn, child_cn, segment_length)

                    child_cn_scores[child_cn] = transition_score + child.cn_score[child_cn]

                min_score = np.amin(child_cn_scores)

                child.cn_backtrack[cn] = np.where(min_score == child_cn_scores)[0]

                node.cn_score[cn] += min_score


def calc_score_recursive_vect(node, leaf_cns, segment_length, cn_max, num_site):
    if node.is_leaf:
        node.cn_score = np.full((num_site, cn_max + 1), np.inf)

        #extract the position of each copy from leaf_cns
        #print("node", node.name)
        #print("list", leaf_cns[node.name])
        #print("score", node.cn_score)

        for site, cp_num in enumerate(leaf_cns[node.name]):
            node.cn_score[site, cp_num] = 0.
        #print("score_after", node.cn_score)

    else:
        node.cn_score = np.zeros((num_site, cn_max + 1))
        for child in node.children:
            calc_score_recursive_vect(child, leaf_cns, segment_length, cn_max, num_site)
            child.cn_backtrack = [[]] * (cn_max + 1)

        for cn in range(0, cn_max + 1):
            #print("--------cn = %d---------"%cn)
            for child in node.children:
                child_cn_scores = np.zeros((num_site, cn_max + 1))
                for child_cn in range(0, cn_max + 1):
                    #print("-------Child_cn = %d-------" %child_cn)
                    transition_score = calculate_transition_score(cn, child_cn, segment_length)
                    #print("Trans score", transition_score)
                    #print("cn_s", child_cn_scores)
                    #print("plus", transition_score + child.cn_score[:, child_cn])
                    child_cn_scores[:, child_cn] = transition_score + child.cn_score[:, child_cn]
                    #print("cn_s_after", child_cn_scores)
                min_score = np.amin(child_cn_scores, axis=1)
                #print("min_score", min_score)
                #child.cn_backtrack[cn] = np.where(min_score == child_cn_scores)[0]

                node.cn_score[:, cn] += min_score


'''
def calc_score_dp(node, leaf_cns, segment_length, cn_max):
    for child in node.leaves:
        #initialize each:
    for n in nodes #BFS reverse order:
        n.cn_score[i] = sum(min(all_child_current))
'''



def calculate_transition_score(parent_cn, child_cn, segment_length):
    if parent_cn == child_cn:
        return 0

    elif parent_cn == 0:
        return np.inf

    else:
        return segment_length * abs(parent_cn - child_cn)


def add_events(tree):
    for child in tree.children:
        child.event = child.cn - tree.cn

        add_events(child)


def annotate_deletion(node, cn):
    for child in node.children:

        child.deletion = False

        for child_cn in child.cn_backtrack[cn]:
            if child_cn < cn:
                child.deletion = True

            annotate_deletion(child, child_cn)


def annotate_copy_number(pos, seg, columns=['major', 'minor']):
    results = []

    for site in seg['sample_id'].unique():
        site_pos = pos[pos['sample_id'] == site]

        site_seg = seg[seg['sample_id'] == site]

        print
        seg.head()
        print
        site_pos.head()
        for chrom in seg.chromosome.unique():
            _pos = site_pos[site_pos['chrom'] == chrom]

            _seg = site_seg[site_seg['chromosome'] == chrom]

            results.append(find_overlapping_segments(_pos, _seg, columns))

    return pd.concat(results)


def find_overlapping_segments(pos, seg, columns):
    seg = seg.sort(['start', 'end'])

    beg_idx = np.searchsorted(seg['start'].values, pos['coord'].values) - 1

    end_idx = np.searchsorted(seg['end'].values, pos['coord'].values)

    mask = (beg_idx == end_idx)

    results = pos.copy()

    for col in columns:
        results[col] = np.nan

        results.loc[mask, col] = seg[col].iloc[end_idx[mask]].values

    return results


def tabulate_cnv_losses(cnvs, nodes, tree):
    # Extract snv chromosome coord ref alt from event id
    nodes['chrom'] = nodes['event_id'].apply(lambda a: a.split(':')[0])
    nodes['coord'] = nodes['event_id'].apply(lambda a: a.split(':')[1]).astype(int)
    nodes['ref'] = nodes['event_id'].apply(lambda a: a.split(':')[2])
    nodes['alt'] = nodes['event_id'].apply(lambda a: a.split(':')[3])

    # Add leaf name to nodes
    node_to_sample_id = dict([(leaf.label, leaf.name) for leaf in tree.leaves])
    node_to_sample_id = pd.Series(node_to_sample_id).reset_index()
    node_to_sample_id.columns = 'node', 'sample_id'
    nodes = nodes.merge(node_to_sample_id, on='node', how='left')

    cn_cols = ['major', 'minor', 'major_is_allele_a']
    nodes = annotate_copy_number(nodes, cnvs, columns=cn_cols)

    nodes['allele_a'] = np.where(nodes['major_is_allele_a'].values == 1, nodes['major'], nodes['minor'])

    nodes['allele_b'] = np.where(nodes['major_is_allele_a'].values == 0, nodes['major'], nodes['minor'])

    # CNV tree with normal node as root of the snv tree
    cnv_tree = dollo.trees.TreeNode()
    cnv_tree.children.append(tree)
    cnv_tree.label = None

    # Create a nodes of possible deletions under minimum evolution
    deletions = list()

    for (chrom, coord, ref, alt), snv in nodes.groupby(['chrom', 'coord', 'ref', 'alt']):
        for allele in ('a', 'b'):
            if snv['allele_' + allele].max() > 10 or snv['allele_' + allele].isnull().any():
                continue

            calc_score_recursive(cnv_tree, snv.set_index('sample_id')['allele_' + allele], 1, 10)

            annotate_deletion(cnv_tree, 1)

            deletion = pd.DataFrame(tree.get_attr_map('deletion').items(), columns=['node', 'deletion'])

            deletion['chrom'] = chrom
            deletion['coord'] = coord
            deletion['ref'] = ref
            deletion['alt'] = alt
            deletion['allele'] = allele

            deletions.append(deletion)

    deletions = pd.concat(deletions, ignore_index=True)

    return deletions