import requests
import logging
import math
from pathlib import Path
import gzip
import hashlib
import tarfile
import re

from tqdm import tqdm

from autometa.config.utilities import DEFAULT_FPATH


# Set up logger
logger = logging.getLogger(__name__)

# --------------------- MD5 Checksum Calculation ---------------------
def calculate_md5(filepath, chunk_size=1024 * 1024):
    """Calculates the MD5 checksum of a file"""
    md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        while chunk := f.read(chunk_size):
            md5.update(chunk)
    return md5.hexdigest()


# --------------------- GTDB Version Handling ---------------------
def get_latest_gtdb_version(host):
    """Fetches the latest GTDB version number from the GTDB server"""
    try:
        response = requests.get(f"https://{host}/releases/latest/VERSION.txt")
        response.raise_for_status()
        version = response.text.splitlines()[0]
        version = version[1:] if version.startswith("v") else version
        return version
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"Failed to fetch GTDB version: {e}")


# --------------------- Taxdump Release URL Fetch ---------------------
def get_gtdb_taxdump_release_url(gtdb_version):
    """Finds the download URL for the GTDB taxdump file for a specific GTDB release"""
    if gtdb_version == "latest":
        raise ValueError(
            "Latest version not supported. Please specify a version number."
        )
    try:
        releases_url = "https://api.github.com/repos/shenwei356/gtdb-taxdump/releases"
        response = requests.get(releases_url)
        response.raise_for_status()
        releases = response.json()
        version_used = str(math.floor(float(gtdb_version)))

        for release in releases:
            if f"r{version_used}" in release["name"]:
                for asset in release["assets"]:
                    if "gtdb-taxdump.tar.gz" in asset["name"]:
                        download_url = asset["browser_download_url"]
                        logger.info(f"Download URL found: {download_url}")
                        return download_url
        logger.error(f"Version R{gtdb_version} not found.")
        return None
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to fetch releases: {e}")
        return None


# --------------------- GTDB Taxdump Download ---------------------
def download_gtdb_taxdump(gtdb_version, outpath, force=False):
    """Downloads the GTDB taxdump file for a specific GTDB release"""
    if not force and Path(outpath).exists():
        logger.info(f"File already exists: {outpath}")
        return outpath
    try:
        download_url = get_gtdb_taxdump_release_url(gtdb_version)
        if download_url:
            response = requests.get(download_url, stream=True)
            response.raise_for_status()
            total_size = int(response.headers.get("content-length", 0))
            chunk_size = 1024
            with tqdm(
                total=total_size, unit="B", unit_scale=True, desc=str(outpath)
            ) as pbar:
                with open(outpath, "wb") as f:
                    for chunk in response.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
            logger.info(f"Download complete. File saved to: {outpath}")
        else:
            logger.error("Download URL was not found.")
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to download GTDB taxdump: {e}")
    except IOError as e:
        logger.error(f"File write error: {e}")
    return outpath


def unpack_gtdb_taxdump(tar_file, gtdb_version, outdir=None, force=False):
    """Extracts the GTDB taxdump file and renames the directory to include the GTDB version"""
    if not outdir:
        outdir = Path(tar_file).parent
    target_dir = f"gtdb-taxdump/R{gtdb_version}"
    new_dir_prefix = f"gtdb_taxdump-version-{gtdb_version}"
    if not force and Path(outdir, new_dir_prefix).exists():
        logger.info(
            f"Directory already exists: {outdir}/{new_dir_prefix}, use --force to overwrite."
        )
        return Path(outdir, new_dir_prefix)
    with tarfile.open(tar_file, "r:gz") as tar:
        members = []
        for member in tar.getmembers():
            if member.name.startswith(target_dir):
                # Adjust the path to rename the folder on extraction
                member.name = member.name.replace(target_dir, new_dir_prefix, 1)
                members.append(member)
        if members:
            tar.extractall(outdir, members=members)
            print(f"Extracted and renamed {target_dir} to {new_dir_prefix} in {outdir}")
        else:
            print(f"Directory {target_dir} not found in the archive.")
    return Path(outdir, new_dir_prefix)


# --------------------- Proteins AA Reps Download with MD5 Verification ---------------------
def download_proteins_aa_reps(host, version, subversion, outpath, force=False):
    """Downloads the GTDB proteins_aa_reps tarball for a specific GTDB release, with MD5 checksum verification"""
    if not force and Path(outpath).exists():
        logger.info(f"File already exists: {outpath}")
        return
    if version == "latest":
        try:
            version = get_latest_gtdb_version(host)
        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to fetch GTDB version number: {e}")
            raise
    logger.info(f"Downloading gtdb_proteins_aa_reps.tar.gz, version {version}")
    try:
        md5sum_url = f"https://{host}/releases/release{version}/{version}.{subversion}/MD5SUM.txt"
        response = requests.get(md5sum_url)
        response.raise_for_status()
        md5sum_lines = response.text.splitlines()
        expected_md5 = None
        filename = f"genomic_files_reps/gtdb_proteins_aa_reps_r{version}.tar.gz"
        for line in md5sum_lines:
            if filename in line:
                expected_md5 = line.split()[0]
                break
        if not expected_md5:
            logger.error(
                f"MD5 checksum for version {version} not found in {md5sum_url}."
            )
            return
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to fetch MD5SUM.txt: {e}")
        return
    url = f"https://{host}/releases/release{version}/{version}.{subversion}/genomic_files_reps/gtdb_proteins_aa_reps_r{version}.tar.gz"
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            logger.info(f"Downloading from {url}")
            total_size = int(r.headers.get("content-length", 0))
            chunk_size = 1024 * 1024
            md5 = hashlib.md5()

            with tqdm(
                total=total_size, unit="B", unit_scale=True, desc=str(outpath)
            ) as pbar:
                with open(outpath, "wb") as f:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
                            md5.update(chunk)
            calculated_md5 = md5.hexdigest()
            if calculated_md5 == expected_md5:
                logger.info(f"MD5 checksum verification passed for {outpath}.")
            else:
                logger.error(
                    f"MD5 checksum verification failed for {outpath}. Expected {expected_md5}, got {calculated_md5}."
                )
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to download the file: {e}")
    except IOError as e:
        logger.error(f"File write error: {e}")
    return outpath


# --------------------- Combined GTDB FASTA Creation ---------------------
def create_combined_gtdb_fasta(tar_file: str, outpath: str, force=False):
    """
    Generate a combined faa file to create the GTDB-t database.

    Parameters
    ----------
    tar_file : str
        The downloaded gtdb_proteins_aa_reps_*.tar.gz file.
    outpath : str
        Path to the combined FASTA file to be written.
    force : bool, optional
        If True, overwrite the output file if it already exists, by default False

    Returns
    -------
    str
        Path to combined faa file. This can be used to make a diamond database.
    """
    # Check if the output file already exists
    if not force and Path(outpath).exists():
        logger.info(f"File already exists: {outpath}")
        return outpath
    # Open the combined output file
    with gzip.open(outpath, "wt") as f_out:
        # Open the tar.gz archive
        with tarfile.open(tar_file, "r:gz") as tar:
            # Initialize tqdm progress bars
            with tqdm(
                desc="Files read", unit="file", position=0, leave=True
            ) as file_pbar, tqdm(
                desc="Sequences written", unit="seq", position=1, leave=True
            ) as seq_pbar:
                # Iterate over members in the tar file
                for member in tar:
                    # Check if the member is a file and ends with .faa.gz
                    if member.isfile() and member.name.endswith(".faa.gz"):
                        # Search for genome accession in the file name
                        genome_acc_search = re.search(
                            r"(GCA_\d+\.\d+|GCF_\d+\.\d+)", member.name
                        )
                        if genome_acc_search:
                            genome_acc = genome_acc_search.group()
                        else:
                            raise ValueError(
                                f"Could not find genome accession for {member.name}"
                            )
                        # Extract and read the content of the .faa.gz file
                        with tar.extractfile(member) as f_in:
                            with gzip.GzipFile(fileobj=f_in) as gz_in:
                                seq_count = (
                                    0  # Initialize sequence counter for the file
                                )
                                for line in gz_in:
                                    line = line.decode("utf-8")
                                    if line.startswith(">"):
                                        seqheader = line.lstrip(">").strip()
                                        outline = f">{genome_acc} {seqheader}\n"
                                        seq_pbar.update(seq_count)
                                    else:
                                        outline = line
                                    f_out.write(outline)
                                    seq_count += 1  # Increment sequence count
                        file_pbar.update(
                            1
                        )  # Update file progress bar after processing each file
    logger.debug(f"Combined GTDB faa file written to {outpath}")
    return outpath


def download_and_format(gtdb_host, gtdb_version, single_dir, force=False):
    """
    Download and format GTDB and NCBI files.

    Parameters
    ----------
    gtdb_host : str
        The GTDB host to download files from.
    gtdb_version : str
        The GTDB version to download.
    single_dir : str
        The single directory to download and format files into.
    dryrun : bool, optional
        If True, only print what would be done, by default False
    """
    if gtdb_version == "latest":
        gtdb_version = get_latest_gtdb_version(gtdb_host)
        logger.info(f"Using 'latest' GTDB version: {gtdb_version}")

    if "." in gtdb_version:
        gtdb_version = gtdb_version.split(".")[0]
        gtdb_subversion = gtdb_version.split(".")[1]
    else:
        gtdb_subversion = "0"
    gtdb_taxdmp_path = Path(single_dir, f"gtdb-taxdump-version-{gtdb_version}.tar.gz")
    # have to rename because GTDB file doesn't have subversion in the name
    aa_reps_path = Path(
        single_dir,
        f"gtdb_proteins_aa_reps-version-{gtdb_version}.{gtdb_subversion}.tar.gz",
    )
    gtdb_taxdmp_path = download_gtdb_taxdump(
        gtdb_version=gtdb_version, outpath=gtdb_taxdmp_path, force=force
    )
    taxdmp_dir = unpack_gtdb_taxdump(
        tar_file=gtdb_taxdmp_path, gtdb_version=gtdb_version, force=force
    )
    aa_reps_path = download_proteins_aa_reps(
        host=gtdb_host,
        version=gtdb_version,
        subversion=gtdb_subversion,
        outpath=aa_reps_path,
        force=force,
    )
    combined_gtdb_fasta = create_combined_gtdb_fasta(
        tar_file=aa_reps_path,
        outpath=Path(
            single_dir,
            f"autometa_formatted_gtdb-version-{gtdb_version}.{gtdb_subversion}.faa.gz",
        ),
        force=force,
    )
    return {
        "gtdb_taxdmp_path": gtdb_taxdmp_path,
        "taxdmp_dir": taxdmp_dir,
        "aa_reps_path": aa_reps_path,
        "combined_gtdb_fasta": combined_gtdb_fasta,
    }



def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="Download GTDB files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version",
        help="GTDB version to download, 'latest' to get the latest version, otherwise specify a version number.",
        default="220",
    )
    parser.add_argument(
        "--host",
        help="GTDB host to download files from.",
        default="data.gtdb.ecogenomic.org",
    )
    parser.add_argument(
        "--outdir",
        help="Directory to save the downloaded files.",
        required=True        
    )
    args = parser.parse_args()
    download_and_format(gtdb_host=args.host, gtdb_version=args.version, single_dir=args.outdir)

if __name__ == "__main__":
    main()
