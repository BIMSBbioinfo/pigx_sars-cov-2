import os
import sys
from snakemake.logging import logger
import download_tarball as dltar

with open(snakemake.log[0], "w") as log_file:
    sys.stdout = log_file
    sys.stderr = sys.stdout

    try:
        logger.info(
            f"Creating database dir at {snakemake.output}...")

        vep_dir = snakemake.output[0]
        os.makedirs(vep_dir)

        dl_url = snakemake.params["dl_url"]

        tar_out = dltar.download_tarball(dl_url)

        tar_out.extractall(vep_dir)

        tar_out.close()

        logger.info(
            f"Download successfull, downloaded vep db files to "
            f"{vep_dir}.")

        logger.info("Done with vep_db.")

    except Exception as e:
        logger.error(e)
        sys.exit(1)
