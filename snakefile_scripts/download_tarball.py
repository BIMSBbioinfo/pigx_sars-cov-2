import re
import requests
import tarfile
import sys
import urllib.request as urlrequest
from snakemake.logging import logger


def download_tarball(dl_url):
    """
    Download .tar.gz archive from url and return a tarfile object to it.
    """

    if not re.search(".tar.gz$", dl_url):
        raise Exception(f"dl_url ({dl_url}) does not point to a '.tar.gz' file.")
    
    if re.match("ftp://*", dl_url):
        logger.info(
            f"Downloading database archive from ftp server at "
            f"{dl_url} to temp file...")

        try:
            dl_resp = urlrequest.urlopen(dl_url)

        except Exception as e:
            logger.error(e)
            sys.exit(1)

        tar_out = tarfile.open(fileobj = dl_resp, mode="r:gz")

    elif re.match("http[s]?://*", dl_url) or not re.match("[a-z]*://", dl_url):
        # FIXME Change http download to also use urllib.request for consistency.
        if not re.match("[a-z]*://", dl_url):
            logger.info(
                f"No protocol identifier in {dl_url}, assuming http/https...")
            dl_url = "http://" + dl_url

        logger.info(
            f"Downloading database archive from {dl_url}...")

        dl_resp = requests.get(dl_url, stream=True)

        if not dl_resp.status_code == 200:
            raise Exception(
                f"Download not successfull, status code {dl_resp.status_code}")

        tar_out = tarfile.open(fileobj=dl_resp.raw, mode="r:gz")

    else:
        raise Exception(
            f"Cannot resolve url protocol in {dl_url}")

    return(tar_out)
