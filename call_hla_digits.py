#!/usr/bin/python
from __future__ import print_function
import sys
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vbseqfile', required=True)
    parser.add_argument('-a', '--allelelist', required=True, help='IMGT HLA Allelelist')
    parser.add_argument('-r', '--mean_rlen', required=True, type=float, help='mean single read length')
    parser.add_argument('-d', '--digit', default=4, type=int, choices=[4 ,6 ,8])
    parser.add_argument('--ispaired', action='store_true', help='if set, twice the mean rlen for depth calculation')
    args = parser.parse_args()
    depth_threshold = 5.0

    alist = pd.read_table(args.allelelist, sep=' ', names=['id', 'subtype'])
    alist['fasta_id'] = 'HLA:' + alist['id']
    alist['gene'] = alist.subtype.map(lambda x: x.split('*')[0])
    alist['subtype_8d'] = alist.subtype.map(lambda x: ':'.join(x.split(':')[:4]))
    alist['subtype_6d'] = alist.subtype.map(lambda x: ':'.join(x.split(':')[:3]))
    alist['subtype_4d'] = alist.subtype.map(lambda x: ':'.join(x.split(':')[:2]))
    alist['subtype_2d'] = alist.subtype.map(lambda x: ':'.join(x.split(':')[:1]))

    default_cn = 2 # fixed

    mean_rlen = args.mean_rlen
    if args.ispaired:
        mean_rlen *= 2

    """
    ID      LENGTH  Z       FPKM    THETA
    HLA:HLA00001    1047    0.00    0.0     0.000000e+00
    HLA:HLA00002    1047    0.00    0.0     0.000000e+00
    HLA:HLA00398    1047    0.00    0.0     0.000000e+00
    HLA:HLA00644    1047    0.00    0.0     0.000000e+00
    HLA:HLA00003    1047    45.74   7877.1781230    1.456456e-06
    """
    tab = pd.read_table(args.vbseqfile)
    tab['mean_depth'] = 1. * mean_rlen * tab['Z'] / tab['LENGTH']
    tab = tab.merge(alist, left_on='ID', right_on='fasta_id')
    del tab['fasta_id']
    del tab['id']

    gene_depths = {}
    for gene, gtab in tab.groupby('gene'):
        gene_depths[gene] = gtab.mean_depth.sum()

    keys = '8d 6d 4d 2d'.split(' ')
    depths = {key: {} for key in keys}

    for key in keys:
        for gkey, gtab in tab.groupby('subtype_' + key):
            depths[key][gkey] = gtab.mean_depth.sum()

    tab['gene_cn'] = tab.gene.map(lambda x: default_cn)
    tab['gene_depth'] = tab.gene.map(lambda x: gene_depths[x])
    tab['ratio_in_gene'] = tab.mean_depth / tab.gene_depth
    for key in keys:
        d = tab['subtype_' + key].map(lambda x: depths[key][x])
        tab['ratio_' + key] = d / tab.gene_depth

    # add genotype
    def add_gt(tab, id_key='ID', ratio_key='ratio_in_gene', gt_key='gt', sel_key='sel'):
        def get_2copy(gtab):
            gtab = gtab.drop_duplicates(id_key, keep='first')
            depths = gtab[ratio_key] * gtab['gene_depth']
            depth_sorted = depths.sort_values(ascending=False)
            depth_tops = depth_sorted[depth_sorted >= depth_threshold]  # filter
            tops = gtab.loc[depth_tops.index]
            if len(tops) >= 2:
                if depth_tops.iloc[0] >= depth_tops.iloc[1] * 2:
                    return ('hom', [tops.iloc[0][id_key]])
                else:
                    return ('het', [tops.iloc[0][id_key], tops.iloc[1][id_key]])
            elif len(tops) == 1:
                if depth_tops.iloc[0] >= depth_threshold * 2:
                    return ('hom', [tops.iloc[0][id_key]])
                else:
                    return ('het', [tops.iloc[0][id_key]])  # partly missing
            else:
                return ('unknown', [])

        gene_gts = {}
        selected_ids = set()
        for gene, gtab in tab.groupby('gene'):
            gt, ids = get_2copy(gtab)
            gene_gts[gene] = gt
            selected_ids.update(set(ids))

        tab[gt_key] = tab.gene.map(lambda x: gene_gts[x])
        tab[sel_key] = tab[id_key].map(lambda x: int(x in selected_ids))

    add_gt(tab)
    for key in keys:
        add_gt(tab, id_key='subtype_' + key, ratio_key='ratio_' + key, gt_key='gt_' + key, sel_key='sel_' + key)

    cond_show = (tab['mean_depth'] > 0)
    tab = tab[cond_show]
    tab = tab.sort_values('gene')
    tab = convert_hla_gt_format(tab)

    geneDic = {}
    for row in tab.query('sel == 1 and digit == @args.digit').itertuples():
        gene, subtype, allele_cn = row.gene, row.subtype, str(row.allele_cn)
        subcn = subtype + '_' + allele_cn
        if gene in geneDic:
            geneDic[gene] = geneDic[gene] + '+' + subcn
        else:
            geneDic[gene] = subcn

    genes = ["A", "B", "C", "E", "F", "G", "DQA1", "DQB1",  "DPB1", "DRB1","L", "MICA", "MICB", "TAP1", "TAP2", "DPA1", "DMA", "DMB", "DOA", "DOB", "DRA", "V",]

    print("Gene\tAllele1\tAllele2")
    for gene in genes:
        sub1, sub2 ='',''
        if gene in geneDic.keys():
            v = geneDic[gene]
            if '+' in v:
                sub1, sub2 = getSub2(v)
            else:
                sub1, allele_cn=v.split('_')
                if allele_cn=='2':
                    sub2 = sub1
        if sub1:
            print(gene, sub1, sub2, sep='\t')

def getSub2(v):
    x1,x2 = sorted(v.split('+'))
    return x1.split('_')[0], x2.split('_')[0]


def convert_hla_gt_format(org_tab):
    digit = 0
    tabs = []
    cols = 'gene digit subtype gene_cn allele_cn allele sel gene_depth ratio refnames'.split(' ')
    digits = [2, 4, 6, 8, 0]

    for gene, tab1 in org_tab.groupby('gene'):
        for digit in digits:
            if digit == 0:
                id_key = 'ID'
                ratio_key = 'ratio_in_gene'
                gt_key = 'gt'
                sel_key = 'sel'
            else:
                id_key = 'subtype_{0}d'.format(digit)
                ratio_key = 'ratio_{0}d'.format(digit)
                gt_key = 'gt_{0}d'.format(digit)
                sel_key = 'sel_{0}d'.format(digit)

            rtab = tab1.groupby(id_key).agg({'ID': lambda x: ','.join(x)}).rename(columns={'ID': 'refnames'}).reset_index()
            tab2 = tab1.drop_duplicates(id_key, keep='first')
            tab2 = tab2.merge(rtab)
            tab2 = tab2.sort_values(ratio_key, ascending=False).reset_index()
            recs = {col: [] for col in cols}
            for an, rec in enumerate(tab2.itertuples(), 1):
                recs['gene'].append(gene)
                recs['digit'].append(digit)
                recs['gene_cn'].append(rec.gene_cn)
                recs['gene_depth'].append(rec.gene_depth)
                recs['allele'].append(an)
                gt = getattr(rec, gt_key)
                if rec.gene_cn == 0:
                    acn = 0
                elif rec.gene_cn == 1 and gt == 'hom':
                    acn = 1
                elif rec.gene_cn == 2 and gt == 'hom':
                    acn = 2
                elif rec.gene_cn == 2 and gt == 'het':
                    acn = 1
                else:
                    acn = 0
                sel = getattr(rec, sel_key)

                recs['allele_cn'].append(acn * sel)
                recs['subtype'].append(getattr(rec, id_key))
                recs['sel'].append(sel)
                recs['ratio'].append(getattr(rec, ratio_key))
                recs['refnames'].append(rec.refnames)
            tab_out = pd.DataFrame(recs, columns=cols)
            tabs.append(tab_out)

    tab = pd.concat(tabs)[cols]
    return tab

if __name__ == '__main__':
    main()

