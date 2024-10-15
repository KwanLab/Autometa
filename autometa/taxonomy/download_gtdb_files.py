import requests
import logging
import math
import os
import glob
logger = logging.getLogger(__name__)



def get_release_url(gtdb_version):
    try:
        # Fetch the list of releases
        releases_url = "https://api.github.com/repos/shenwei356/gtdb-taxdump/releases"
        response = requests.get(releases_url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        releases = response.json()
        # Search through each release and its assets for a matching file (gtdb-taxdump-R220)
        for release in releases:
            for asset in release['assets']:
                # Make version a lower integer (since decimal isn't present in gtdb-taxdump)
                version_used = str(math.floor(float(gtdb_version)))
                if f"R{version_used}" in asset['name']:
                    download_url = asset['browser_download_url']
                    logger.info(f"Download URL found: {download_url}")
                    return download_url
        logger.error(f"Version R{gtdb_version} not found.")
        return None
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to fetch releases: {e}")
        return None


def download_gtdb_taxdump(gtdb_version, outdir, force=False):
    if not force:
        # Check if the file already exists
        filepath = os.path.join(outdir, f"gtdb-taxdump-R{gtdb_version}.tar.gz")
        if os.path.exists(filepath):
            logger.info(f"File already exists: {filepath}")
            return
    try:
        download_url = get_release_url(gtdb_version)
        if download_url:
            # Download the GTDB taxdump file
            response = requests.get(download_url)
            response.raise_for_status()  # Raise an exception for HTTP errors
            filepath = os.path.join(outdir, f"gtdb-taxdump-R{gtdb_version}.tar.gz")
            with open(filepath, 'wb') as f:
                f.write(response.content)
            logger.info(f"Download complete. File saved to: {filepath}")
        else:
            logger.error("Download URL was not found.")
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to download GTDB taxdump: {e}")
    except IOError as e:
        logger.error(f"File write error: {e}")

def download_proteins_aa_reps(host,version,outdir,force=False):
    if not force:
        # Stop if the file already exists, use glob to find the file because version not known
        filepath = os.path.join(outdir, "gtdb_proteins_aa_reps_r*.tar.gz")
        if len(glob.glob(filepath)) > 0:
            logger.info(f"File already exists: {filepath}")
            return
    # if version is integer string:
    if str(version).isdigit():
        version = str(version)
    elif version == "latest":
        # don't download from latest url, find the latest version number then download that version
        try:
            response = requests.get(f"https://{host}/releases/latest/VERSION.txt")
            response.raise_for_status()
            version = response.text.splitlines()[0]
            # version = "v220" -> "220"
            version = version.removeprefix("v")
        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to fetch GTDB version number: {e}")
            raise
    logger.info(f"Downloading gtdb_proteins_aa_reps.tar.gz, version {version}")
    url = f"https://{host}/releases/release{version}/{version}.0/genomic_files_reps/gtdb_proteins_aa_reps_r{version}.tar.gz"
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        logger.info(f"Downloading from {url}")
        with open(os.path.join(outdir, f"gtdb_proteins_aa_reps_r{version}.tar.gz"), 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)


