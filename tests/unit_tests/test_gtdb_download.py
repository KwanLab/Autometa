import pytest
from unittest import mock
import requests
import os
from io import BytesIO
import logging
from autometa.taxonomy.download_gtdb_files import get_latest_gtdb_version, get_gtdb_taxdump_release_url, download_gtdb_taxdump, download_proteins_aa_reps

# Mock logger to avoid unnecessary logging during tests
logger = logging.getLogger(__name__)


@pytest.fixture
def mock_requests_get():
    with mock.patch('requests.get') as mock_get:
        yield mock_get


def test_get_latest_gtdb_version_success(mock_requests_get):
    # Mock a successful request response
    mock_response = mock.Mock()
    mock_response.status_code = 200
    mock_response.text = "v220\n"
    mock_requests_get.return_value = mock_response

    host = "data.gtdb.ecogenomic.org"
    version = get_latest_gtdb_version(host)
    
    assert version == "220"
    mock_requests_get.assert_called_once_with(f"https://{host}/releases/latest/VERSION.txt")


def test_get_latest_gtdb_version_fail(mock_requests_get):
    # Mock a failed request
    mock_requests_get.side_effect = requests.exceptions.RequestException("Error occurred")
    
    host = "data.gtdb.ecogenomic.org"
    with pytest.raises(requests.exceptions.RequestException):
        get_latest_gtdb_version(host)
    

def test_get_gtdb_taxdump_release_url_success(mock_requests_get):
    # Mock a response from GitHub API with a valid release
    mock_response = mock.Mock()
    mock_response.status_code = 200
    mock_response.json.return_value = [
        {
            'name': 'r220',
            'assets': [{'name': 'gtdb-taxdump.tar.gz', 'browser_download_url': 'https://example.com/download'}]
        }
    ]
    mock_requests_get.return_value = mock_response

    gtdb_version = "220"
    url = get_gtdb_taxdump_release_url(gtdb_version)
    
    assert url == 'https://example.com/download'
    mock_requests_get.assert_called_once_with("https://api.github.com/repos/shenwei356/gtdb-taxdump/releases")


def test_get_gtdb_taxdump_release_url_not_found(mock_requests_get):
    # Mock a response with no matching version
    mock_response = mock.Mock()
    mock_response.status_code = 200
    mock_response.json.return_value = [
        {'name': 'r219', 'assets': [{'name': 'gtdb-taxdump.tar.gz'}]}
    ]
    mock_requests_get.return_value = mock_response

    gtdb_version = "220"
    url = get_gtdb_taxdump_release_url(gtdb_version)
    
    assert url is None
    mock_requests_get.assert_called_once_with("https://api.github.com/repos/shenwei356/gtdb-taxdump/releases")


def test_download_gtdb_taxdump_file_exists():
    # Mock the file already existing
    with mock.patch('os.path.exists', return_value=True):
        with mock.patch('requests.get') as mock_get:
            download_gtdb_taxdump("220", "/some/dir", force=False)
            mock_get.assert_not_called()

def test_download_gtdb_taxdump_success(mock_requests_get):
    # Mock a successful file download
    mock_response = mock.Mock()
    mock_response.status_code = 200
    mock_response.headers = {'content-length': '1024'}
    mock_response.iter_content = mock.Mock(return_value=[b'chunk1', b'chunk2'])
    # Mock the JSON response for the releases
    mock_response.json.return_value = [
        {
            'name': 'r220',
            'assets': [{'name': 'gtdb-taxdump.tar.gz', 'browser_download_url': 'https://example.com/download'}]
        }
    ]
    mock_requests_get.return_value = mock_response
    with mock.patch('os.path.exists', return_value=False):
        with mock.patch('builtins.open', mock.mock_open()) as mock_file:
            with mock.patch('autometa.taxonomy.download_gtdb_files.tqdm') as mock_tqdm:
                download_gtdb_taxdump("220", "/some/dir", force=False)
                mock_requests_get.assert_called()
                mock_file.assert_called_once_with('/some/dir/gtdb-taxdump-R220.tar.gz', 'wb')
                mock_tqdm.assert_called_once()



def test_download_proteins_aa_reps_success(mock_requests_get):
    # Mock successful file download for proteins_aa_reps
    mock_response = mock.Mock()
    mock_response.status_code = 200
    mock_response.headers = {'content-length': '1024'}
    mock_response.iter_content = mock.Mock(return_value=[b'chunk1', b'chunk2'])
    # Add context manager support
    mock_requests_get.return_value.__enter__.return_value = mock_response
    mock_requests_get.return_value.__exit__ = mock.Mock()
    with mock.patch('os.path.exists', return_value=False):
        with mock.patch('builtins.open', mock.mock_open()) as mock_file:
            # Mock tqdm in the specific module where it's used (update the import path accordingly)
            with mock.patch('autometa.taxonomy.download_gtdb_files.tqdm') as mock_tqdm:
                download_proteins_aa_reps("data.gtdb.ecogenomic.org", "220", "/some/dir", force=False)
                mock_requests_get.assert_called()
                mock_file.assert_called_once_with('/some/dir/gtdb_proteins_aa_reps-R220.tar.gz', 'wb')
                mock_tqdm.assert_called_once()
