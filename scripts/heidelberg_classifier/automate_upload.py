import requests
from utils.log import get_console_logger, get_file_logger
from utils.output import unique_output_dir
import os
from settings import HEIDELBERG_CLASSIFIER_CONFIG, DATA_DIR
import pandas as pd
from scripts.heidelberg_classifier import api
from time import sleep

logger = get_console_logger(__name__)


def batch_upload(
    samplesheet_fn,
    idat_dir,
    include=None,
    exclude=None,
    array_type='EPIC',
    sample_type='FFPE DNA',
    include_subdir=True,
    n_retry=3,
    wait_between_retries=5
):
    """

    :param samplesheet_fn:
    :param idat_dir:
    :param include:
    :param exclude:
    :param array_type: One of 'EPIC', '450K'
    :param sample_type: One of 'FFPE DNA', 'KRYO DNA'
    :param include_subdir: If True (default), look for the pair of idat files in a subdirectory within idat_dir, named
    by the sample_id. This is the format for the data returned by the UCL facility, but not for GEO downloads.
    :return:
    """
    analysis_params = {
        'chipType': array_type,
        'sampleType': sample_type
    }
    obj = api.Heidelberg()
    samples = api.read_samplesheet(samplesheet_fn)
    if include is not None:
        samples = samples.loc[samples.Sample_Name.isin(include)]
    if exclude is not None:
        samples = samples.loc[~samples.Sample_Name.isin(exclude)]

    logger.info("Uploading %d samples", len(samples))
    for i, row in samples.iterrows():
        n = 1
        sample_id = str(row.Sentrix_ID)  # may be an int
        name = "%s;%s;%s" % (base_dir, sample_id, row.Sample_Name)
        sample_pos = row.Sentrix_Position
        if include_subdir:
            file = os.path.join(idat_dir, sample_id, "%s_%s_Red.idat" % (sample_id, sample_pos))
        else:
            file = os.path.join(idat_dir, "%s_%s_Red.idat" % (sample_id, sample_pos))
        success = False
        while success is False:
            try:
                # if this is a retry, log back in first
                if n > 1:
                    obj.establish_session()
                obj.submit_sample(name, file, **analysis_params)
                success = True
                logger.info("Success: %s", name)
            except Exception as exc:
                logger.exception("Failed to process sample %s", name)
                if n < n_retry:
                    n += 1
                    logger.info("Retrying (attempt %d) in %d seconds", n, wait_between_retries)
                    logger.info("Retrying %s. Attempt %d.", name, n)
                    sleep(wait_between_retries)

        if not success:
            logger.error("Upload failed: %s", name)
            # raise AttributeError("Failed to upload sample %s after %d retries" % (name, n_retry))

    return obj


"""
Maintained for purely historic reasons: this is how I classified the TCGA raw data.
The files to upload are found within the subdirs. It is assumed that each subdir has only 2 idats
"""
# the_login_url = 'https://www.molecularneuropathology.org/mnp/authenticate'
# the_submission_url = 'https://www.molecularneuropathology.org/mnp/sample/add'
# indir = '/home/gabriel/python_outputs/gdc-nih_paired_methylation_gene_expr.2'
# cids = [
#     # '9b9654c9-5a64-4a8e-b1f0-987437e6012e',
#     # '4644a15f-3115-4f27-86b2-b092419431e1',
#     'c2f7f72f-13ae-4efc-9cda-69ca037b4498',
#     '28da11ec-f46a-4abc-b4c8-bbb784419c71',
#     '9d1d9ae5-4640-4060-a207-c9c5b0e905e1',
#     '872abc8a-6c1f-4114-b993-7d0327fb38bd',
#     # '3a3fc890-1985-4353-861b-dc3abfb364b1',
#     '8da3103e-3e6c-4176-a583-d5fe5e60601e',
#     '7a650a2c-bc3f-4e0c-820e-4f492977107b',
#     '74ce7e7f-e214-4b6e-8551-114a153ab7fb',
#     'af055b98-be6a-4012-ac64-f1b6539b23d0',
#     '51c3409d-5df7-4720-83db-0e6de065f82c',
#     'c29d73c0-c885-4105-bf74-38e9178e71c9',
#     'd73ca945-b60b-4b5c-ac53-7ef2602a4951',
#     '153c2442-ea61-4b5e-8c5c-c71d287c6055',
#     '7338e476-10ff-444a-a16f-429d28355f65',
#     '3f960d3b-a58c-43d0-a8a4-f3555b399c9d',
#     'e876dd29-68b0-4bf1-83d1-488c40068a35',
#     '34f216fb-09dc-48af-9e09-7d12cc07c1f7',
#     '20bad0d5-3135-49a4-a119-fa7a1e56fd1b',
#     '752a6e21-9a91-417f-beba-bfdf331d5cac',
#     '8e6b9705-e9d2-4fbc-ac03-35c1a5115eef',
#     'd6de8d1c-e5ff-45cd-b53e-f943d2578713',
#     # '2c5fa2d4-8e35-42e4-8bca-9fb3371a19c8',
#     'c04657d2-b71b-4402-82fa-02747cce331d',
#     'a28b256e-e085-467d-bdba-5c39718012d7',
#     'dfa3ef71-7006-4c6d-81c4-e0c0c7b75c26',
#     'd558203d-3f70-42d4-accb-008ace280f48',
#     'e6e2d4c3-d37c-4de3-ac68-b301800770f0',
#     '8df94420-c736-4a36-b001-5dbe3dff17f0',
#     '69d56f2d-6924-409b-9d1e-c8d69b400270',
#     'a262928c-e20a-4c02-8114-1227e05c43e1',
#     '9348e446-0e43-4210-b07e-c534980cbf62',
#     '74139255-a635-4c87-814d-3dd04ed630a8',
#     'd59ffef5-f39b-4627-a7c7-6920b99c2408',
#     'd1547b99-3c96-4c62-8261-5111bcf860a9',
#     '5c984433-33cf-42fc-b3ba-511efcdcab19',
#     'c184c3ca-7ad3-4202-b108-cb9fd5f5d947',
#     '46bfebf4-0ef6-4348-81dc-d7d3cb52c08f',
#     '93ed7a2b-b0cb-4a84-871f-5c34a0b6a640',
#     '5327e899-a20d-4571-8236-98454bad574e',
#     'f4fe7c02-be19-4929-8638-960e5776494c',
#     '4c42dc4e-66b7-40bf-9fbe-b92543248198',
#     '8d2e88d9-d8d0-4c42-8aa2-205a788dea58',
#     '884f867b-4a8b-4b67-8fe4-ab3f068be84e',
#     'de7b7cac-f094-4d59-8651-e991e34ea093',
#     'd9fcf27f-48e3-43fd-a8cf-3db246cf4221',
#     '6ccd57e5-9215-46fa-b3be-e408a0f424dc',
#     'c970240f-03e2-4395-85be-5eb81a69b710',
#     'c9e4f0bf-48fc-4600-9e03-314bc575273f'
# ]
#
# field_defaults = {
#     'age': '',
#     'location': '',
#     'diagnosis': '',
#     'notes': '',
#     'gender': 'NA',
#     'chipType': 'NA',
#     'sampleType': 'KRYO DNA',
# }
#
# s = requests.Session()
# resp = s.post(the_login_url, data=HEIDELBERG_CLASSIFIER_CONFIG)
# if resp.status_code != 200:
#     logger.error("Login failed")
#
# for cid in cids:
#     logger.info("Case ID %s", cid)
#     fulldir = os.path.join(indir, cid)
#     idats = [t for t in os.listdir(fulldir) if 'idat' in t.lower()]
#     if len(idats) != 2:
#         logger.error("Expected 2 idat files in directory %s. Found %d.", fulldir, len(idats))
#         continue
#     payload = {
#         'name': cid,
#         'diagnosis': 'GBM',
#         'chipType': '450k',
#     }
#     files = {
#         'file1': open(os.path.join(fulldir, idats[0]), 'rb'),
#         'file2': open(os.path.join(fulldir, idats[1]), 'rb'),
#     }
#     # set default values for all else
#     for f, v in field_defaults.items():
#         payload.setdefault(f, v)
#     resp = s.post(the_submission_url, data=payload, files=files)
#     if resp.status_code != 200:
#         logger.error("Upload failed for sample %s: %s", cid, resp.content)
#     else:
#         sleep(1)


if __name__ == '__main__':
    """
    Automated upload of samples to the Heidelberg classifier site
    """
    inputs = [
        #{'base_dir': '2016-06-10_brandner',
        # 'chipType': 'EPIC',
        # 'sampleType': 'FFPE DNA'},
        #{'base_dir': '2016-09-21_dutt',
        # 'chipType': 'EPIC',
        # 'sampleType': 'FFPE DNA'},
        # {'base_dir': '2016-12-19_ucl_genomics',
        # 'chipType': 'EPIC',
        # 'sampleType': 'KRYO DNA'},
        # {'base_dir': '2017-01-17_brandner',
        #  'chipType': 'EPIC',
        #  'sampleType': 'FFPE DNA'},
        # {'base_dir': '2017-02-09_brandner',
        #  'chipType': 'EPIC',
        #  'sampleType': 'FFPE DNA'},
        # {'base_dir': '2017-05-12',
        #  'chipType': 'EPIC',
        #  'sampleType': 'KRYO DNA'},
        # {'base_dir': '2017-08-23',
        #  'chipType': 'EPIC',
        #  'sampleType': 'KRYO DNA'},
        # {'base_dir': '2017-09-19',
        #  'chipType': 'EPIC',
        #  'sampleType': 'KRYO DNA'},
        {'base_dir': 'GSE92462_450K',
         'chipType': '450K',
         'sampleType': 'KRYO DNA'},
    ]
    # include_subdir = True
    include_subdir = False
    n_retry = 3
    wait_between_retries = 5  # seconds
    outdir = unique_output_dir('heidelberg_bulk_upload', reuse_empty=True)
    flog = get_file_logger('heidelberg_bulk_upload', os.path.join(outdir, 'automated_upload.log'))
    # include = [
    #     'GBM018 P10 DNA 8/11/2016 CLEANED',
    #     'GBM019 P4 DNA 8/11/2016 CLEANED',
    #     'GBM024 P9 DNA     13/10/2016',
    #     'GBM026 P8 DNA 24/10/2016',
    #     'GBM031 P4 DNA     13/10/2016'
    # ]
    include = [
        'GSC80',
        'GSC164',
        'GSC64',
        'GSC76',
        'GSC102',
        'GSC6',
        'GSC14',
        'GSC10',
        'GSC39',
        'GSC43',
        'GSC108',
        'GSC59',
        'GSC12',
        'GSC84',
        'GSC61',
        'GSC22',
        'GSC187',
        'GSC182',
        'GSC154',
        'GSC184',
        'GSC143',
        'GSC195',
        'H1228',
    ]
    exclude = None

    api_objs = {}
    obj = api.Heidelberg()

    for d in inputs:
        base_dir = d.pop('base_dir') # remaining dictionary values are kwargs for the uploader
        print "Project %s" % base_dir
        sample_fn = os.path.join(DATA_DIR, 'methylation', base_dir, 'samplesheet.csv')
        indir = os.path.join(DATA_DIR, 'methylation', base_dir, 'idat')

        samples = api.read_samplesheet(sample_fn)
        if include is not None:
            samples = samples.loc[samples.Sample_Name.isin(include)]
        if exclude is not None:
            samples = samples.loc[~samples.Sample_Name.isin(exclude)]

        print "Uploading %d samples" % len(samples)
        flog.info("Project %s. %d samples.", base_dir, len(samples))
        for i, row in samples.iterrows():
            n = 1
            # sample_id = str(int(row.Sentrix_ID))  # may be an int
            sample_id = str(row.Sentrix_ID)  # may be an int
            name = "%s;%s;%s" % (base_dir, sample_id, row.Sample_Name)
            sample_pos = row.Sentrix_Position
            if include_subdir:
                file = os.path.join(indir, sample_id, "%s_%s_Red.idat" % (sample_id, sample_pos))
            else:
                file = os.path.join(indir, "%s_%s_Red.idat" % (sample_id, sample_pos))
            success = False
            while success is False:
                try:
                    # if this is a retry, log back in first
                    if n > 1:
                        obj.establish_session()
                    obj.submit_sample(name, file, **d)
                    success = True
                    print "Uploaded %s" % name
                    flog.info("Success: %s", name)
                except Exception as exc:
                    print "Failed to process sample %s: %s" % (name, repr(exc))
                    if n < n_retry:
                        n += 1
                        print "Retrying (attempt %d) in %d seconds" % (n, wait_between_retries)
                        flog.info("Retrying %s. Attempt %d.", name, n)
                        sleep(wait_between_retries)

            if not success:
                flog.error("Upload failed: %s", name)
                # raise AttributeError("Failed to upload sample %s after %d retries" % (name, n_retry))

        api_objs[base_dir] = obj.submitted
        print "Done"
