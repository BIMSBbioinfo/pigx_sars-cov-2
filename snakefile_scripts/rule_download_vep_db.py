import os
import sys
import requests
import urllib.request as urlrequest
import tarfile
import re
from snakemake.logging import logger

with open(snakemake.log[0], "w") as log_file:
    sys.stdout = log_file
    sys.stderr = sys.stdout

    try:
        logger.info(
            f"Creating database dir at {snakemake.output}...")

        vep_dir = snakemake.output[0]
        os.makedirs(vep_dir)

        dl_url = snakemake.params["dl_url"]

        temp_dl_file = os.path.join(vep_dir, "archive.temp")

        if re.match("ftp://*", dl_url):
            logger.info(
                f"Downloading database archive from ftp server at "
                f"{dl_url} to temp file...")

            urlrequest.urlretrieve(
                dl_url, temp_dl_file)

            tar_out = tarfile.open(temp_dl_file)

        elif re.match("http[s]?://*", dl_url):
            logger.info(
                f"Downloading database archive from {dl_url}...")

            dl_resp = requests.get(dl_url, stream=True)

            if not dl_resp.status_code == 200:
                logger.error(
                    f"Download not successfull, status code {dl_resp.status_code}")
                sys.exit(1)

            tar_out = tarfile.open(fileobj=dl_resp.raw, mode="r:gz")

        tar_out.extractall(vep_dir)

        tar_out.close()

        if os.path.exists(temp_dl_file):
            os.remove(temp_dl_file)

        logger.info(
            f"Download successfull, downloaded vep db files to "
            f"{vep_dir}.")

        logger.info("Done with vep_db.")

    except Exception as e:
        logger.error(e)
        sys.exit(1)
