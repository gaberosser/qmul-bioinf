import pandas as pd
import numpy as np
from scipy import stats
from scipy.cluster import hierarchy
import os
from settings import DATA_DIR
from matplotlib import pyplot as plt
import seaborn as sns
import multiprocessing as mp
from sklearn.decomposition import PCA
from plotting import heatmap, clustering

from utils import setops, output
from .cpt_cmyc import standardise

if __name__ == '__main__':
    myc_probe = 3115514
    myc_gene = 'MYC'
    corr_method = 'pearson'
    cross_corr_threshold = 0.5
    within_corr_threshold = 0.75
    alpha = 0.05

    # pre-determined outlier
    outlier = 'CP62-T'

    # these genes were validated, so we'd like to retain them!
    validated_genes = [
        'CXCL1',
        'CXCL2',
        'IL6',
        'IL1RL2',
    ]

    outdir = output.unique_output_dir("cpt_myc_leverage")

    fn = os.path.join(DATA_DIR, 'exon_array', 'GSE60982', 'raw', 'GSE60892_HuEx-ALL40-Upload-Transposed.txt.gz')
    rma_data = pd.read_csv(fn, sep='\t', comment='#', header=0, index_col=0)

    ann_fn = os.path.join(DATA_DIR, 'exon_array', 'GSE60982', 'HuEx-1_0-st-v2.na36.hg19.probeset.csv.gz')
    ann = pd.read_csv(ann_fn, sep=',', comment='#', header=0, index_col=0, na_values='---')

    # annotation files
    # latest version from Affymetrix / Thermofisher
    # every probe in this annotation is defined as 'core'
    ann = ann.loc[rma_data.index]
    ann = ann.loc[~ann.gene_assignment.isnull()]
    ann.loc[:, 'gene_assignment'] = ann.loc[:, 'gene_assignment'].str.replace(r' /+ ', '|')
    the_symbols = ann.gene_assignment.str.split('|').apply(lambda x: x[1])

    # get a standardised dataset
    the_data = standardise(rma_data)

    # we'll skip computing probes and just declare them in a list (taken from the previous code)
    # this is only for one MYC probe (3115514), but will probably be fine
    keep_probes = [2321858, 2321860, 2324576, 2324774, 2324792, 2329499, 2329518, 2330690, 2330696, 2330697, 2330703,
                   2330704, 2333830, 2334195, 2334202, 2334207, 2334214, 2349021, 2363785, 2363795, 2363797, 2370361,
                   2370362, 2370363, 2370370, 2371182, 2375017, 2378075, 2378076, 2379141, 2379142, 2379147, 2379148,
                   2379149, 2379150, 2379151, 2379152, 2379153, 2379154, 2388560, 2388570, 2390051, 2390057, 2391012,
                   2392962, 2392974, 2392988, 2392990, 2392992, 2392996, 2392998, 2398543, 2402915, 2414960, 2414961,
                   2414962, 2414963, 2415087, 2415089, 2415093, 2415101, 2435075, 2440120, 2443542, 2448384, 2448387,
                   2448395, 2448397, 2451937, 2452207, 2471979, 2474270, 2484372, 2497036, 2497038, 2497046, 2497054,
                   2522257, 2522314, 2530717, 2530718, 2530722, 2532028, 2544791, 2547469, 2553630, 2556219, 2562443,
                   2563664, 2563666, 2563679, 2563682, 2563684, 2563688, 2563703, 2565101, 2565105, 2565106, 2565110,
                   2565120, 2565122, 2565127, 2571487, 2571511, 2571512, 2571517, 2571519, 2571520, 2576078, 2604014,
                   2605736, 2634115, 2634116, 2634118, 2634127, 2634128, 2634130, 2634132, 2643168, 2647745, 2650366,
                   2651997, 2654968, 2654982, 2666961, 2666962, 2669931, 2672197, 2672199, 2672201, 2672741, 2672744,
                   2672753, 2682438, 2691872, 2696048, 2702339, 2702344, 2702348, 2702353, 2702361, 2702367, 2702368,
                   2702369, 2705760, 2705787, 2711672, 2711691, 2715636, 2729854, 2731335, 2731338, 2731339, 2731341,
                   2731354, 2731355, 2731356, 2731383, 2731385, 2731387, 2731388, 2731389, 2768360, 2773389, 2773395,
                   2773396, 2773437, 2773438, 2773439, 2773443, 2773445, 2773598, 2774832, 2776002, 2777588, 2782239,
                   2821257, 2827672, 2827691, 2827692, 2827693, 2864122, 2866707, 2866709, 2875350, 2875351, 2875360,
                   2878079, 2878283, 2887320, 2895252, 2899468, 2900203, 2900215, 2905418, 2909783, 2909784, 2909787,
                   2909789, 2922765, 2924592, 2924594, 2925849, 2927511, 2927513, 2927515, 2927516, 2927518, 2927519,
                   2927520, 2927521, 2927522, 2927523, 2927524, 2927526, 2927527, 2927528, 2927529, 2927530, 2927531,
                   2927532, 2927533, 2931604, 2934115, 2934123, 2948170, 2948171, 2957539, 2958516, 2962858, 2963721,
                   2964717, 2979197, 2979252, 2982090, 2982335, 2992594, 2992596, 2992597, 2992601, 2992602, 2994569,
                   2994706, 3001106, 3012820, 3016149, 3016150, 3016153, 3016154, 3016158, 3016161, 3016163, 3016166,
                   3016168, 3016170, 3016172, 3036930, 3046447, 3061626, 3061627, 3061629, 3061630, 3061632, 3061633,
                   3061634, 3061635, 3063051, 3064374, 3066865, 3082980, 3110057, 3115017, 3115021, 3115022, 3115514,
                   3115515, 3115522, 3115523, 3115524, 3117388, 3127589, 3127590, 3129744, 3140229, 3140230, 3140232,
                   3144255, 3146686, 3147541, 3147542, 3194752, 3194763, 3214455, 3219226, 3219228, 3219230, 3219233,
                   3219235, 3219242, 3219821, 3219856, 3222008, 3225333, 3225406, 3225407, 3225408, 3225410, 3225411,
                   3225416, 3225417, 3225418, 3226201, 3226203, 3227113, 3252045, 3252047, 3252049, 3252052, 3252056,
                   3252058, 3252061, 3252062, 3252064, 3260385, 3262016, 3267325, 3274429, 3292171, 3292267, 3293343,
                   3294960, 3299581, 3299582, 3299975, 3299986, 3302417, 3319946, 3319952, 3319956, 3319965, 3319970,
                   3319971, 3319974, 3319978, 3319979, 3320134, 3320136, 3320137, 3320139, 3322809, 3323083, 3323484,
                   3326683, 3326686, 3326694, 3326700, 3326705, 3326711, 3326712, 3329263, 3334273, 3335521, 3346555,
                   3346563, 3346564, 3346567, 3346574, 3350855, 3353765, 3354888, 3355861, 3355869, 3362266, 3369445,
                   3369491, 3376223, 3376259, 3376261, 3377250, 3377994, 3377996, 3378759, 3378760, 3380975, 3391172,
                   3392050, 3396220, 3405775, 3410657, 3415241, 3442817, 3442819, 3442824, 3442827, 3442866, 3442872,
                   3442874, 3442883, 3443869, 3443870, 3443872, 3443874, 3443880, 3448105, 3449409, 3452666, 3453361,
                   3453683, 3455151, 3455262, 3456612, 3458590, 3461845, 3471290, 3472375, 3475647, 3476756, 3476761,
                   3476768, 3485801, 3519324, 3519327, 3522664, 3533560, 3561040, 3561043, 3561044, 3571554, 3571555,
                   3591666, 3592404, 3592407, 3622210, 3645568, 3645576, 3662116, 3662137, 3662144, 3662145, 3662154,
                   3662193, 3662195, 3662252, 3662253, 3662255, 3676344, 3685241, 3695345, 3695351, 3695353, 3706127,
                   3714747, 3714751, 3718169, 3718170, 3718173, 3718175, 3718176, 3718564, 3718566, 3718932, 3718935,
                   3718937, 3718982, 3718985, 3718986, 3719026, 3719028, 3719030, 3719031, 3719371, 3719525, 3720048,
                   3721688, 3728954, 3728956, 3728957, 3744113, 3744116, 3744963, 3746588, 3754012, 3754017, 3754139,
                   3762585, 3764583, 3764586, 3764588, 3764589, 3764590, 3764594, 3764595, 3764596, 3764598, 3764604,
                   3764607, 3764609, 3764610, 3764611, 3764614, 3764615, 3764617, 3764618, 3764619, 3764620, 3764624,
                   3764625, 3764626, 3764628, 3764629, 3764630, 3764631, 3764633, 3764634, 3764635, 3764636, 3764637,
                   3764640, 3764641, 3764642, 3776524, 3776534, 3776544, 3776545, 3776546, 3802967, 3816512, 3816515,
                   3816517, 3816519, 3816520, 3816521, 3816522, 3816523, 3816524, 3816525, 3816950, 3817562, 3818390,
                   3820444, 3820445, 3820446, 3820450, 3821895, 3821897, 3822221, 3823532, 3832904, 3832979, 3832986,
                   3836284, 3836285, 3836286, 3836288, 3836289, 3836290, 3838005, 3838006, 3838007, 3838262, 3841264,
                   3842337, 3851450, 3852797, 3860112, 3861906, 3864558, 3864559, 3864560, 3864563, 3864564, 3864577,
                   3864829, 3865838, 3871267, 3890662, 3890663, 3895625, 3910261, 3912099, 3914782, 3919024, 3923761,
                   3934920, 3945237, 3945238, 3945241, 3945243, 3947031, 3955649, 3957167, 3957168, 3957170, 3957381,
                   3957387, 3957389, 3969486, 3969490, 3972974, 3974884, 3975543, 3975817, 3975818, 3975822, 3978175,
                   3978176, 3978179, 3978181, 3978184, 3978187, 3978189, 3986090, 3989762, 3989780, 3990738, 3992165,
                   3992416, 3995305, 4000513, 4000514, 4000521, 4000522, 4000531, 4002159, 4002161, 4007646, 4007740,
                   4009114, 4018735, 4018739, 4018746, 4021277, 4024109, 4026577, 4027604, 4027609, 4040065, 4040067,
                   4040070, 4040071]

    # now we want to show how the outlier affects the correlation measure
    cxcl1_probes = [2731383, 2731385, 2731387, 2731388, 2731389]
    cxcl2_probes = [2773437, 2773438, 2773439, 2773443, 2773445]
    il6_probes = [2992594, 2992596, 2992597, 2992601, 2992602]
    il1rl2_probes = [2497036, 2497038, 2497046, 2497054]

    genes = {
        'CXCL1': cxcl1_probes,
        'CXCL2': cxcl2_probes,
        'IL6': il6_probes,
        'IL1RL2': il1rl2_probes,
    }

    myc_probe_data = the_data.loc[myc_probe]
    myc_probe_data_censor = myc_probe_data.drop(outlier)
    myc_probe_data_outlier = myc_probe_data.loc[[outlier]]
    x_lr = np.array([np.floor(myc_probe_data.min()), np.ceil(myc_probe_data.max())])

    # now we'll generate scatterplots for 4 probes in each of these cases
    for nm, arr in genes.items():
        fig, axs = plt.subplots(2, 2, sharex=True, figsize=(10.5, 7.5))
        for i, pid in enumerate(arr[:4]):
            ax = axs.flat[i]
            x_probe_data = the_data.loc[pid]
            fit_all = stats.linregress(myc_probe_data, x_probe_data)
            x_probe_data_censor = x_probe_data.drop(outlier)
            fit_cen = stats.linregress(myc_probe_data_censor, x_probe_data_censor)

            ax.scatter(myc_probe_data_censor.values, x_probe_data_censor.values, c='k')

            y_lr = fit_all.intercept + fit_all.slope * x_lr
            ax.plot(x_lr, y_lr, 'r--', label='Pearson r=%.2f' % fit_all.rvalue)

            ax.scatter(myc_probe_data_outlier.values, x_probe_data.loc[[outlier]].values, c='r')
            y_lr = fit_cen.intercept + fit_cen.slope * x_lr
            ax.plot(x_lr, y_lr, 'k--', label='Pearson r=%.2f' % fit_cen.rvalue)

            ax.set_ylabel('%s probe %d RMA expression' % (nm, pid))
            ax.set_xlabel('MYC probe %d RMA expression' % myc_probe)
            ax.legend(loc='upper left')

        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "probe_level_difference_%s.png" % nm), dpi=200)

    # repeat at the gene level
    # redefine probes now so that we have all of them, not just those correlating with the first MYC probe

    genes = {
        'CXCL1': [2731383, 2731385, 2731387, 2731388, 2731389],
        # 'CXCL2': [2773437, 2773438, 2773439, 2773443, 2773445],
        'IL6': [2992592, 2992594, 2992596, 2992597, 2992601, 2992602],
        'IL1RL2': [2497035, 2497036, 2497038, 2497045, 2497046, 2497047, 2497051, 2497054],
        'DSC3': [3802967],
    }

    myc_gene_data = the_data.loc[[3115514, 3115515, 3115522, 3115523, 3115524]].mean(axis=0)
    myc_gene_data_censor = myc_gene_data.drop(outlier)
    myc_gene_data_outlier = myc_gene_data.loc[[outlier]]

    x_lr = np.array([np.floor(myc_gene_data.min()), np.ceil(myc_gene_data.max())])

    gene_data = dict([
        (nm, the_data.loc[arr].mean(axis=0)) for nm, arr in genes.items()
    ])
    gene_data = pd.DataFrame.from_dict(gene_data).transpose()

    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(10.5, 7.5))
    for i, nm in enumerate(gene_data.index):
        x_gene_data = gene_data.loc[nm]
        ax = axs.flat[i]
        fit_all = stats.linregress(myc_gene_data, x_gene_data)
        x_gene_data_censor = x_gene_data.drop(outlier)
        fit_cen = stats.linregress(myc_gene_data_censor, x_gene_data_censor)

        ax.scatter(myc_gene_data_censor.values, x_gene_data_censor.values, c='k')

        y_lr = fit_all.intercept + fit_all.slope * x_lr
        ax.plot(x_lr, y_lr, 'r--', label='Pearson r=%.2f' % fit_all.rvalue)

        ax.scatter(myc_gene_data_outlier.values, x_gene_data.loc[[outlier]].values, c='r')
        y_lr = fit_cen.intercept + fit_cen.slope * x_lr
        ax.plot(x_lr, y_lr, 'k--', label='Pearson r=%.2f' % fit_cen.rvalue)

        ax.set_ylabel('%s RMA expression' % nm)
        ax.set_xlabel('MYC RMA expression')
        if nm == 'DSC3':
            ax.legend(loc='lower left')
        else:
            ax.legend(loc='upper left')

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "gene_level_difference.png"), dpi=200)