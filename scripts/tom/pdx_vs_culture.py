from methylation import loader
from scripts.hgic_final import consts
import hgic_consts


if __name__ == "__main__":
    """
    Here we assess the similarity of GIC cultures (in vitro) at both early and late passages and those established
    after xenografting (ex vivo).
    We also have bulk tumour samples from human patients (FFPE) and mouse xenografts (FFPE PDX).

    This is all with DNA methylation data (EPIC array).

    Except for FFPE, this is for one line (only): 019. May consider including the others to demonstrate how different
    lines/bulk vary.
    """
    norm_method = 'swan'
    pdx_bulk_samples = [
        'SM18_108A GBM19Luc PDX1',
        'SM18_119A GBM19Luc PDX2'
    ]
    gic_late_samples = [
        'GBM19Luc P12',
        'GBM19Luc PDX1 P3',
        'GBM19Luc PDX2 P2',
    ]

    # load all relevant data
    our_gic_obj = loader.load_by_patient(consts.PIDS, include_control=False, samples=consts.S1_METHYL_SAMPLES_GIC, norm_method=norm_method)
    our_ffpe_obj = loader.load_by_patient(consts.PIDS, type='ffpe', include_control=False, norm_method=norm_method)
    pdx_bulk = loader.load_reference('2018-12-14', norm_method=norm_method, samples=pdx_bulk_samples)
    gic_late = loader.load_reference('2018-12-06', norm_method=norm_method, samples=gic_late_samples)

    # add patient ID to our FFPE
    our_ffpe_obj.meta.insert(0, 'patient_id', [hgic_consts.NH_ID_TO_PATIENT_ID_MAP[t] for t in our_ffpe_obj.meta.index])


