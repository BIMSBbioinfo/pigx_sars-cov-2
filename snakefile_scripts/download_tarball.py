import re
import requests
import tempfile
import tarfile
import urllib.request as urlrequest


def download_tarball(dl_url):
    """
    Download .tar.gz archive from url and return a tarfile object to it.
    """

    if not re.search(".tar.gz$", dl_url):
        Exception(f"dl_url ({dl_url}) does not point to a '.tar.gz' file.")
    
    if re.match("ftp://*", dl_url):
        print(
            f"Downloading database archive from ftp server at "
            f"{dl_url} to temp file...")

        temp_dl_file = tempfile.TemporaryFile()

        urlrequest.urlretrieve(
            dl_url, temp_dl_file)

        tar_out = tarfile.open(temp_dl_file)

    elif re.match("http[s]?://*", dl_url):
        print(
            f"Downloading database archive from {dl_url}...")

        dl_resp = requests.get(dl_url, stream=True)

        if not dl_resp.status_code == 200:
            Exception(
                f"Download not successfull, status code {dl_resp.status_code}")

        tar_out = tarfile.open(fileobj=dl_resp.raw, mode="r:gz")

    else:
        Exception(
            f"Cannot resolve url protocol in {dl_url}")

    return(tar_out)
