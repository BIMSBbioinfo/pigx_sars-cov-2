import os
import sys
import subprocess
import requests
import tarfile
from snakemake.logging import logger

with open(snakemake.log[0], "w") as log_file:
    sys.stdout = log_file
    sys.stderr = sys.stdout

    try:
        krona_dir = snakemake.output[0]

        if snakemake.params["use_prebuilt"]:
            logger.info(
                f"Creating database dir at {krona_dir}...")

            os.makedirs(krona_dir)

            dl_url = snakemake.params["dl_url"]

            logger.info(
                f"Downloading prebuilt database archive from {dl_url}...")

            dl_resp = requests.get(dl_url, stream=True)

            if not dl_resp.status_code == 200:
                logger.error(
                    f"Download not successfull, status code {dl_resp.status_code}")
                sys.exit(1)

            tar_out = tarfile.open(fileobj=dl_resp.raw, mode="r:gz")

            tar_out.extractall(krona_dir)

            tar_out.close()

            logger.info(
                f"Download successfull, downloaded krona db files to "
                f"{krona_dir}.")

        else:
            logger.info(
                f"Letting krona download its current version into "
                f"{krona_dir}")

            krona_update_tax_script = snakemake.params["krona_update_tax_script"]

            krona_update_tax_call = [
                krona_update_tax_script, krona_dir]

            logger.info(krona_update_tax_call)

            subprocess.run(
                krona_update_tax_call,
                stdout=log_file,
                stderr=log_file
            )

        logger.info("Done with krona_db.")

    except Exception as e:
        logger.error(e)
        sys.exit(1)
