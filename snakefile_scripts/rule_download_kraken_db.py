import os
import sys
import shutil
import requests
import tarfile
import subprocess
from snakemake.logging import logger

with open(snakemake.log[0], "w") as log_file:
    sys.stdout = log_file
    sys.stderr = sys.stdout

    try:
        logger.info(
            f"Creating database dir at {snakemake.output}...")

        kraken_dir = snakemake.output[0]
        os.makedirs(kraken_dir)

        db_url = snakemake.params["dl_url"]

        logger.info(
            f"Downloading database archive from {db_url}...")

        dl_resp = requests.get(db_url, stream=True)

        if not dl_resp.status_code == 200:
            logger.error(
                f"Download not successfull, status code {dl_resp.status_code}")
            sys.exit(1)

        tar_out = tarfile.open(fileobj=dl_resp.raw, mode="r:gz")

        downsample_db = snakemake.params["downsample_db"]
        max_db_size = snakemake.params["max_db_size"]

        if downsample_db:
            logger.info(
                f"Downsampling database to {max_db_size}B. "
                f"This requires additional work.")

            for member in tar_out:
                if member.name != "hash.k2d":
                    tar_out.extract(member, path=kraken_dir)

            tar_out.close()

            taxonomy_dl_call = [
                snakemake.params["kraken2_build"],
                "--use-ftp",
                "--download-taxonomy",
                "--skip-map",
                "--db", kraken_dir]

            logger.info(
                "Downloading taxonomy..."
            )

            logger.info(taxonomy_dl_call)

            subprocess.run(
                taxonomy_dl_call,
                stdout=log_file,
                stderr=log_file
            )

            # Manually create library dummy dir. Somehow kraken2-build still
            # expects that.
            os.makedirs(os.path.join(kraken_dir, "library"))

            db_build_call = [
                *snakemake.params["kraken2_build"].split(" "),
                "--build",
                "--db", kraken_dir,
                "--max-db-size", str(max_db_size)]

            logger.info(
                "Building database with size-limit...")

            logger.info(db_build_call)

            subprocess.run(
                db_build_call,
                stdout=log_file,
                stderr=log_file
            )

            # After building, the taxonomy subdir only takes up space.
            shutil.rmtree(os.path.join(kraken_dir, "taxonomy"))

        else:
            logger.info(
                "No downsampling requested, taking downloaded dir as is.")

            tar_out.extractall(path=kraken_dir)

            logger.info(
                f"Download successfull, downloaded kraken2 files to "
                f"{kraken_dir}.")

        tar_out.close()

        logger.info("Done with kraken_db.")

    except Exception as e:
        logger.error(e)
        sys.exit(1)
