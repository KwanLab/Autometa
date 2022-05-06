import sys
from autometa import __version__, console_scripts

import argparse

citation = """
APA:

Miller, I. J., Rees, E. R., Ross, J., Miller, I., Baxa, J., Lopera, J., Kerby, R. L.,
Rey, F. E., & Kwan, J. C. (2019). Autometa: automated extraction of microbial genomes
from individual shotgun metagenomes. Nucleic Acids Research, 47(10).
https://doi.org/10.1093/nar/gkz148

BibTeX:

@article{
    Miller_Autometa_automated_extraction_2019,
    author = {Miller, Ian J. and Rees, Evan R. and Ross, Jennifer and Miller, Izaak and Baxa, Jared and Lopera, Juan and Kerby, Robert L. and Rey, Federico E. and Kwan, Jason C.},
    doi = {10.1093/nar/gkz148},
    journal = {Nucleic Acids Research},
    number = {10},
    title = {{Autometa: automated extraction of microbial genomes from individual shotgun metagenomes}},
    url = {https://github.com/KwanLab/Autometa},
    volume = {47},
    year = {2019}
}
"""


def main():
    parser = argparse.ArgumentParser(
        description="Describe Autometa citation & version."
        "No arguments will list the available autometa commands, docs and code information"
    )
    parser.add_argument(
        "-V", "--version", help="Print autometa version", action="store_true"
    )
    parser.add_argument(
        "-C",
        "--citation",
        help="Print autometa citation (APA and BibTex)",
        action="store_true",
    )
    args = parser.parse_args()

    if args.version:
        print(f"autometa: {__version__}")
        sys.exit(0)

    if args.citation:
        print(citation)
        sys.exit(0)

    print("Autometa Commands\n\t│")
    commands_header = "\t├──> "
    commands_body = "\n\t├──> ".join(list(console_scripts)[:-1])
    commands_footer = f"└──> {list(console_scripts)[-1]}"

    print(f"{commands_header}{commands_body}\n\t{commands_footer}")
    print(
        "\nRun 'autometa --version' or 'autometa --citation' for respective info"
        "\nDocs: https://autometa.readthedocs.io"
        "\nCode: https://github.com/KwanLab/Autometa"
    )


if __name__ == "__main__":
    main()
