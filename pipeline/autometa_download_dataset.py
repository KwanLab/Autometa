# in response to issue #110
# Intended: fetch data similarly to scikit-learn API
# pulling data from google drive folder with simulated or synthetic communities

# use gdown to download data from google drive to output directory specified by the user
# create a dictionary of the databases in the google drive
# allow the user to call them based on size (eg '78', '156'...)
# allow the user to specify <some/directory>
# find that corresponding file and download it to <some/directory>

# goal: autometa-download-dataset --community 78 --output <some/directory>

# prepare dependencies
import gdown
import argparse

# take in commands that user input
# including test file for now
parser = argparse.ArgumentParser(prog='autometa-download-dataset', description='Download a simulated community file from google drive to a specified directory')
parser.add_argument('--community', 
                    help='specify a size of simulated community in MB', 
                    choices=['78', '156', '312', '625', '1250', '2500', '5000', '10000', 'test'], 
                    required=True)
parser.add_argument('--output', 
                    help='specify the directory to download the file', 
                    required=True)
args = parser.parse_args()

# provide list of database options as a dictionary with file_ids from google
simulated = {
    'test': '1fy3M7RnS_HGSQVKidCy-rAwXuxldyOOv',
    '78': '15CB8rmQaHTGy7gWtZedfBJkrwr51bb2y',
    '156': '13bkwFBIUhdWVWlAmVCimDODWF-7tRxgI',
    '312': '1qyAu-m6NCNuVlDFFC10waOD28j15yfV-',
    '625': '1FgMXSD50ggu0UJbZd1PM_AvLt-E7gJix',
    '1250': '1KoxwxBAYcz8Xz9H2v17N9CHOZ-WXWS5m',
    '2500': '1wKZytjC4zjTuhHdNUyAT6wVbuDDIwk2m',
    '5000': '1IX6vLfBptPxhL44dLa6jePs-GRw2XJ3S',
    '10000': '1ON2vxEWC5FHyyPqlfZ0znMgnQ1fTirqG'
}

# construct file id into a url to put into gdown
file_id = simulated[args.community]
url = f'https://drive.google.com/uc?id={file_id}'

# download the specified file with gdown
gdown.download(url, args.output)