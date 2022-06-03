#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2022 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
Shaurya Chanana, Izaak Miller, Jason C. Kwan

This file is part of Autometa.

Autometa is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Autometa is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Autometa. If not, see <http://www.gnu.org/licenses/>.
COPYRIGHT

Template Script for Autometa Modules

Template Description:
Concise sentence describing functionality of script

Template Format:
0. - Shebang python env definition
1. - Two lines following comment block script description
2. - One line separating different import types
3. - Two lines separating imports and logger
4. - Two lines separating logger and algorithm functions
5. - Main function
6. - if __name__ == '__main__' clause
7. - Argparse
8. - Logging aliased to logger in 5.
"""


import logging
from typing import Literal, Union
import pandas as pd


# Two lines to separate imports and beginning of script
logger = logging.getLogger(__name__)
# Two lines to separate function definitions and logger


def format_taxon_binning(
    df: pd.DataFrame, sample_id: str, version: str = "0.9.0"
) -> str:
    """Format autometa taxon binning results to be compatible with specified BioBox `version`.

    Parameters
    ----------

    taxon_binning : Union[str,pd.DataFrame]
        Path to (to-be-formatted) taxon_binning results
        or `pd.DataFrame(index_col=range(0, n_rows), columns=[contig, taxid, ...])`

    sample_id : str
        Sample identifier, not the generating user or program name.
        It MUST match the regular expression `[A-Za-z0-9\._]+.`

    version : str, optional
        Biobox version to format results, by default `"0.9"`

    Returns
    -------
    str
        formatted results ready to be written to a file path

    Raises
    -------
    NotImplementedError
        Specified `version` is not implemented
    TypeError
        `taxon_binning` must be a path to a taxon results file or pandas
    ValueError
        `taxon_binning` does not contain the required columns
    """
    biobox_version = f"@VERSION:{version}"
    biobox_sampleid = f"@SAMPLEID:{sample_id}"
    if "taxid" not in df.columns:
        raise ValueError(
            f"taxon_binning results require columns 'contig' and 'taxid'. contains: {df.columns}"
        )
    outcols = ["@@SEQUENCEID", "TAXID"]
    if "cluster" in df.columns:
        df = df.rename(columns={"cluster": "BINID"})
        outcols.append("BINID")
    df = df.rename(columns={"contig": "@@SEQUENCEID", "taxid": "TAXID"})[outcols]
    df_str = df.to_csv(sep="\t", index=False, header=True)
    return f"{biobox_version}\n{biobox_sampleid}\n{df_str}"


def format_genome_binning(
    df: pd.DataFrame, sample_id: str, version: str = "0.9.0"
) -> str:
    """Format autometa genome binning results to be compatible with specified BioBox `version`.

    Parameters
    ----------

    genome_binning : Union[str,pd.DataFrame]
        Path to (to-be-formatted) genome_binning results
        or `pd.DataFrame(index_col=range(0, n_rows), columns=[contig, taxid, ...])`

    sample_id : str
        Sample identifier, not the generating user or program name.
        It MUST match the regular expression `[A-Za-z0-9\._]+.`

    version : str, optional
        Biobox version to format results, by default `"0.9"`

    Returns
    -------
    str
        formatted results ready to be written to a file path

    Raises
    -------
    NotImplementedError
        Specified `version` is not implemented
    TypeError
        `genome_binning` must be a path to a taxon results file or pandas
    ValueError
        `genome_binning` does not contain the required columns
    """
    biobox_version = f"@VERSION:{version}"
    biobox_sampleid = f"@SAMPLEID:{sample_id}"
    if "cluster" not in df.columns:
        raise ValueError(
            f"genome_binning results require columns 'contig' and 'cluster'. contains: {df.columns}"
        )
    outcols = ["@@SEQUENCEID", "BINID"]
    if "taxid" in df.columns:
        df = df.rename(columns={"taxid": "TAXID"})
        outcols.append("TAXID")
    df = df.rename(columns={"contig": "@@SEQUENCEID", "cluster": "BINID"})[outcols]
    df_str = df.to_csv(sep="\t", index=False, header=True)
    return f"{biobox_version}\n{biobox_sampleid}\n{df_str}"


def format_profiling(df: pd.DataFrame, sample_id: str, version: str) -> str:
    raise NotImplementedError(
        "Autometa currently does not perform profiling, are you sure this is the format you need?"
    )


def get_biobox_format(
    predictions: Union[str, pd.DataFrame],
    sample_id: str,
    results_type: Literal["profiling", "genome_binning", "taxon_binning"],
    version: str,
) -> str:
    versions = {"0.9.0"}
    if version not in versions:
        raise NotImplementedError(
            f"{version} not implemented. Available versions: {', '.join(versions)}"
        )

    if isinstance(predictions, str):
        df = pd.read_csv(predictions, sep="\t")
    elif isinstance(predictions, pd.DataFrame):
        df = predictions.copy()
    else:
        raise TypeError(
            f"predictions must be a path to a file or a DataFrame object! Given: {type(predictions)}"
        )

    formatter_dispatcher = {
        "profiling": format_profiling,
        "genome_binning": format_genome_binning,
        "taxon_binning": format_taxon_binning,
    }
    if results_type not in formatter_dispatcher:
        raise ValueError(
            f"{results_type} not in formatters! {', '.join(formatter_dispatcher.keys())}"
        )

    logger.info(
        f"Formatting predictions to biobox {results_type} format (version:{version})"
    )

    formatter = formatter_dispatcher[results_type]
    return formatter(df, sample_id=sample_id, version=version)


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="""
        Format Autometa results to biobox format for compatibility with CAMI.

        Note: All results tables must contain a 'contig' column and either 'taxid', or 'cluster' column.

        bioboxes formatted columns will be written for both if both are within the provided results.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--sample-predictions",
        help="Path of autometa table containing relevant `results-type` columns",
        required=True,
    )
    parser.add_argument(
        "--sample-id",
        help="CAMI Sample ID corresponding to `sample-predictions`",
        required=True,
    )
    parser.add_argument(
        "--results-type",
        help="Type of results for formatter to convert",
        choices=["profiling", "genome_binning", "taxon_binning"],
        required=True,
    )
    parser.add_argument(
        "--bioboxes-version",
        help="bioboxes binning output format. For more info see: https://github.com/bioboxes/rfc/blob/4bb19a633a6a969c2332f1f298852114c5f89b1b/data-format/binning.mkd",
        default="0.9.0",
    )
    parser.add_argument(
        "--output", help="Path to write biobox formatted results", required=True
    )
    args = parser.parse_args()

    biobox_results = get_biobox_format(
        predictions=args.sample_predictions,
        sample_id=args.sample_id,
        results_type=args.results_type,
        version=args.bioboxes_version,
    )

    with open(args.output, "w") as fh:
        fh.write(biobox_results)

    logger.info(f"Wrote bioboxes formatted results to {args.output}")


if __name__ == "__main__":
    main()
