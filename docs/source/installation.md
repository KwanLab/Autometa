(installation-page)=

# Installation

Currently Autometa package installation is supported by [conda] and [docker].
For installation using conda, we suggest downloading [miniconda].

:::{attention}
If you are only trying to run the Autometa workflow, you should start at {ref}`Getting Started` before proceeding.
:::

## Direct installation (Quickest)

1. Install [miniconda]

2. Create a new environment with `autometa` installed:

   > ```bash
   > conda create -c bioconda -n autometa autometa
   > ```

3. Activate `autometa` environment:

   > ```
   > conda activate autometa
   > ```

## Install from source (using make)

Download and install [miniconda]. Now run the following commands:

```bash
# Navigate to the directory where you would like to clone Autometa
cd $HOME

# Clone the Autometa repository
git clone https://github.com/KwanLab/Autometa.git

# Navigate into the cloned repository
cd Autometa

# create autometa conda environment
make create_environment

# activate autometa conda environment
conda activate autometa

# install autometa source code in autometa environment
make install
```

:::{note}
You can see a list of all available make commands by running `make` without any other arguments.
:::

## Install from source (full commands)

Download and install [miniconda]. Now run the following commands:

```bash
# Navigate to the directory where you would like to clone Autometa
cd $HOME

# Clone the Autometa repository
git clone https://github.com/KwanLab/Autometa.git

# Navigate into the cloned repository
cd Autometa

# Construct the autometa environment from autometa-env.yml
conda env create --file=autometa-env.yml

# Activate environment
conda activate autometa

# Install the autometa code base from source
python -m pip install . --ignore-installed --no-deps -vv
```

## Building the Docker image

You can build a docker image for your clone of the Autometa repository.

1. Install [Docker]
2. Run the following commands

```bash
# Navigate to the directory where you need to clone Autometa
cd $HOME

# Clone the Autometa repository
git clone https://github.com/KwanLab/Autometa.git

# Navigate into the cloned repository
cd Autometa

# This will tag the image as jasonkwan/autometa:<your current branch>
make image

# (or the full command from within the Autometa repo)
docker build . -t jasonkwan/autometa:`git branch --show-current`
```

## Testing Autometa

You can also check the installation using autometa's built-in unit tests.
This is not at all necessary and is primarily meant for development and debugging purposes.
To run the tests, however, you'll first need to install the following packages and download the test dataset.

```bash
# Activate your autometa conda environment
conda activate autometa

# List all make options
make

# Install dependencies for test environment
make test_environment

# Download test_data.json for unit testing to tests/data/
make unit_test_data_download
```

You can now run different unit tests using the following commands:

```bash
# Run all unit tests
make unit_test

# Run unit tests marked with entrypoint
make unit_test_entrypoints

# Run unit tests marked with WIP
make unit_test_wip
```

:::{note}
As a shortcut you can also create the test environment and run **all** the unit tests using `make unit_test` command.
:::

For more information about the above commands see the {ref}`Contributing Guidelines` page.
Additional unit tests are provided in the test directory. These are designed to aid in future development of autometa.

[anaconda]: https://www.anaconda.com/
[conda]: https://docs.conda.io/en/latest/
[docker]: https://www.docker.com/
[miniconda]: https://docs.conda.io/en/latest/miniconda.html