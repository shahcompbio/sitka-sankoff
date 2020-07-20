import matplotlib
from matplotlib.colors import ListedColormap
from scgenome import refgenome
import pandas as pd
import numpy as np
import seaborn

color_reference = {0: '#3182BD', 1: '#9ECAE1', 2: '#CCCCCC', 3: '#FDCC8A', 4: '#FC8D59', 5: '#E34A33', 6: '#B30000',
                   7: '#980043', 8: '#DD1C77', 9: '#DF65B0', 10: '#C994C7', 11: '#D4B9DA'}


def get_cn_cmap(cn_data):
    min_cn = int(cn_data.min())
    max_cn = int(cn_data.max())
    assert min_cn - cn_data.min() == 0
    assert max_cn - cn_data.max() == 0
    color_list = []
    for cn in range(min_cn, max_cn + 1):
        if cn > max(color_reference.keys()):
            cn = max(color_reference.keys())
        color_list.append(color_reference[cn])
    return ListedColormap(color_list)


def plot_cell_cn_profile(ax, cn_data, value_field_name, cn_field_name=None, max_cn=13, chromosome=None, s=5):
    """ Plot copy number profile on a genome axis
    Args:
        ax: matplotlib axis
        cn_data: copy number table
        value_field_name: column in cn_data to use for the y axis value

    Kwargs:
        cn_field_name: state column to color scatter points
        max_cn: max copy number for y axis
        chromosome: single chromosome plot
        s: size of scatter points
    The cn_data table should have the following columns (in addition to value_field_name and
    optionally cn_field_name):
        - chr
        - start
        - end
    """
    chromosome_info = refgenome.info.chromosome_info[['chr', 'chromosome_start', 'chromosome_end']].copy()
    cat =  [str(i) for i in range (1, 23)]
    cat.append("X")
    cat.append("Y")
    print("cat", cat)
    chromosome_info['chr'] = pd.Categorical(chromosome_info['chr'], categories=cat)
    plot_data = cn_data.merge(chromosome_info)
    plot_data = plot_data[plot_data['chr'].isin(refgenome.info.chromosomes)]
    plot_data['start'] = plot_data['start'] + plot_data['chromosome_start']
    plot_data['end'] = plot_data['end'] + plot_data['chromosome_start']

    if cn_field_name is not None:
        ax.scatter(
            plot_data['start'], plot_data[value_field_name],
            c=plot_data[cn_field_name], s=s,
            cmap=get_cn_cmap(plot_data[cn_field_name].astype(int).values),
        )
    else:
        ax.scatter(
            plot_data['start'], plot_data[value_field_name], s=s,
        )

    if chromosome is not None:
        chromosome_length = refgenome.info.chromosome_info.set_index('chr').loc[chromosome, 'chromosome_length']
        chromosome_start = refgenome.info.chromosome_info.set_index('chr').loc[chromosome, 'chromosome_start']
        chromosome_end = refgenome.info.chromosome_info.set_index('chr').loc[chromosome, 'chromosome_end']
        xticks = np.arange(0, chromosome_length, 1e7)
        xticklabels = ['{0:d}M'.format(int(x / 1e6)) for x in xticks]
        xminorticks = np.arange(0, chromosome_length, 1e6)
        ax.set_xlabel(f'chromosome {chromosome}')
        ax.set_xticks(xticks + chromosome_start)
        ax.set_xticklabels(xticklabels)
        ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(xminorticks + chromosome_start))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlim((chromosome_start, chromosome_end))

    else:
        ax.set_xlim((-0.5, refgenome.info.chromosome_end.max()))
        ax.set_xlabel('chromosome')
        ax.set_xticks([0] + list(refgenome.info.chromosome_end.values))
        ax.set_xticklabels([])
        ax.xaxis.tick_bottom()
        ax.yaxis.tick_left()
        ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(refgenome.info.chromosome_mid))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(refgenome.info.chromosomes))

    if cn_field_name is not None:
        ax.set_ylim((-0.5, max_cn))
        ax.set_yticks(np.arange(0, max_cn, 2))

    if chromosome is not None:
        seaborn.despine(ax=ax, offset=10, trim=False)
    else:
        seaborn.despine(ax=ax, offset=10, trim=True)

    return chromosome_info
