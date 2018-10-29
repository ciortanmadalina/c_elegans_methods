import numpy as np
import pandas as pd
from scipy import sparse, io
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
from collections import Counter
from IPython.display import clear_output, Image, display
from collections import OrderedDict
from scipy import stats
from statsmodels.stats.multitest import multipletests


def ttest_scores(df, clusters, top_genes = 5):
    """
    This method returns arrays of t test scores, p values and top genes
    """
    #Array keeping the cluster names (keys)
    clusterNames= list(zip(*sorted(Counter(clusters).items())))[0]
    # Array keeping the size of clusters
    ns = list(zip(*sorted(Counter(clusters).items())))[1]
    rankings_gene_scores = []
    rankings_gene_names = []
    rankings_gene_logfoldchanges = []
    rankings_gene_pvals = []
    rankings_gene_pvals_adj = []
    for igroup in clusterNames:
        idx = np.where(clusters == igroup)[0]
        idx_rest = np.where(clusters != igroup)[0]

        mean_group = df.iloc[idx].mean(axis = 0)
        var_group = df.iloc[idx].var(axis = 0)

        mean_rest = df.iloc[idx_rest].mean(axis = 0)
        var_rest = df.iloc[idx_rest].var(axis = 0)

        ns_group = len(idx)
        ns_rest = len(idx_rest)

        denominator = np.sqrt(var_group/ns_group + var_rest/ns_rest)
        denominator[np.flatnonzero(denominator == 0)] = np.nan
        scores = (mean_group - mean_rest) / denominator #Welch t-test
        mean_rest[mean_rest == 0] = 1e-9  # set 0s to small value
        foldchanges = (mean_group + 1e-9) / mean_rest
        scores[np.isnan(scores)] = 0
        #https://en.wikipedia.org/wiki/Welch%27s_t-test
        #Get p-values and degree of freedom
        denominator_dof = (np.square(var_group) / (np.square(ns_group)*(ns_group-1))) + (
            (np.square(var_rest) / (np.square(ns_rest) * (ns_rest - 1))))
        denominator_dof[np.flatnonzero(denominator_dof == 0)] = np.nan
        dof = np.square(var_group/ns_group + var_rest/ns_rest) / denominator_dof # dof calculation for Welch t-test
        dof[np.isnan(dof)] = 0
        pvals = stats.t.sf(abs(scores), dof)*2 # *2 because of two-tailed t-test

        pvals[np.isnan(pvals)] = 1  # set Nan values to 1 to properly convert using Benhjamini Hochberg
    #     _, pvals_adj, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')

        scores_sort_idx = np.argsort(scores)[::-1][:top_genes]

        rankings_gene_scores.append(scores[scores_sort_idx])
    #     rankings_gene_logfoldchanges.append(np.log2(np.abs(foldchanges[scores_sort_idx])))
        rankings_gene_names.append(df.columns[scores_sort_idx])
        rankings_gene_pvals.append(pvals[scores_sort_idx])
    #     rankings_gene_pvals_adj.append(pvals_adj[scores_sort_idx])

    #     print(scores[scores_sort_idx], df.columns[scores_sort_idx])
    
    return rankings_gene_scores, rankings_gene_names, rankings_gene_pvals


def plot_ttest(clusters, rankings_gene_scores, rankings_gene_names, rankings_gene_pvals, figname = '1.pdf'):
    #Array keeping the cluster names (keys)
    clusterNames= list(zip(*sorted(Counter(clusters).items())))[0]
    
    n_panels_x = 4
    n_panels_y = np.ceil(len(clusterNames) / n_panels_x).astype(int)
    # n_panels_y = np.ceil(20 / n_panels_x).astype(int)
    print(n_panels_x, n_panels_y)

    fig = plt.figure(figsize = (n_panels_x * 4, n_panels_y * 4))
    left = 0.2/n_panels_x
    bottom = 0.13/n_panels_y
    gs = gridspec.GridSpec(nrows=n_panels_y,
                           ncols=n_panels_x,
                           left=left,
                           right=1-(n_panels_x-1)*left-0.01/n_panels_x,
                           bottom=bottom,
                           top=1-(n_panels_y-1)*bottom-0.1/n_panels_y,
                           wspace=0.22,
                           hspace=0.4)

    ax0 = None
    ymin = np.Inf
    ymax = -np.Inf
    for count, group_name in enumerate(clusterNames):
        ax = fig.add_subplot(gs[count])
        gene_names = rankings_gene_names[count]
        scores = rankings_gene_scores[count]
        pvals =rankings_gene_pvals[count]
        plt.grid()
        plt.scatter(np.arange(len(gene_names)),scores, s=30, alpha = 0.5)
        for ig, g in enumerate(gene_names):
            gene_name = gene_names[ig]
            ax.text(
                ig,scores[ig],
                f"pval {round(pvals[ig],5)}" , verticalalignment='bottom',
                horizontalalignment='center', fontsize=10)
        plt.xticks(np.arange(len(gene_names)), gene_names, rotation = 20)


        ax.set_title('Cluster {} vs. rest'.format(group_name))
        if count >= n_panels_x * (n_panels_y - 1):
            ax.set_xlabel('top ranking genes')

        # print the 'score' label only on the first panel per row.
        if count % n_panels_x == 0:
            ax.set_ylabel('t score')

        ax.set_xlim(-0.9, ig + 1-0.1)
        ymin = np.min(scores)
        ymax = np.max(scores)
        ymax += 0.3*(np.max(scores)-np.min(scores))
        ax.set_ylim(ymin, ymax)
        ns = list(zip(*sorted(Counter(clusters).items())))[1]
        ns = np.array(ns)
        fig.suptitle(
            f'Total clusters {len(ns)}, valid(#>=10) {len(ns[ns>=10])}, invalid {len(ns[ns<10])}, outliers {np.sum(ns[ns<10])}',
            size=12 )
        plt.savefig(f'report/{figname}')
