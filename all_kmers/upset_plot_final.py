#!/usr/bin/env python3
"""
Create an UpSet plot of k-mer presence across pangenome references.

Expected DuckDB database:
    /gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/
    output_files/kmer_presence.db

Expected table:
    kmer_presence

Expected columns:
    kmer
    African_pan
    HPRC2
    CPG
    Arab-PG
    KSPG

Presence columns should contain 0/1 values.

The intersection counts are calculated directly inside DuckDB, so the
script does not load every individual k-mer into pandas.

Compatible with Python 3.8.
"""

# ============================================================
# Install missing packages
# ============================================================

import importlib
import subprocess
import sys


REQUIRED_PACKAGES = {
    "duckdb": "duckdb",
    "pandas": "pandas",
    "matplotlib": "matplotlib",
    "upsetplot": "upsetplot",
}


def install_missing_packages():
    """Install required packages that are not already available."""
    for module_name, package_name in REQUIRED_PACKAGES.items():
        try:
            importlib.import_module(module_name)

        except ImportError:
            print(
                "Package '{}' was not found. Installing...".format(
                    package_name
                )
            )

            subprocess.check_call(
                [
                    sys.executable,
                    "-m",
                    "pip",
                    "install",
                    "--user",
                    package_name,
                ]
            )


install_missing_packages()


# ============================================================
# Imports
# ============================================================

import warnings
from pathlib import Path
from typing import Tuple

import duckdb
import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet


# ============================================================
# Configuration
# ============================================================

DATABASE = Path(
    "/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/"
    "output_files/kmer_presence.db"
)

TABLE_NAME = "kmer_presence"

OUTPUT_DIRECTORY = Path(
    "/gpfs/home/hce24xau/scratch/gen_kmers/data/all_kmers/"
    "output_files"
)

OUTPUT_PNG = OUTPUT_DIRECTORY / "kmer_upset_plot.png"

OUTPUT_PDF = OUTPUT_DIRECTORY / "kmer_upset_plot.pdf"

OUTPUT_COUNTS = (
    OUTPUT_DIRECTORY / "kmer_upset_intersection_counts.csv"
)


# Display label mapped to the DuckDB column name.
#
# The order here controls the initial set order.
SETS = {
    "African Pan": "African_pan",
    "HPRC2": "HPRC2",
    "CPG": "CPG",
    "Arab-PG": "Arab-PG",
    "KSPG": "KSPG",
}


# Minimum exact intersection size to retain.
#
# Use None to retain all observed intersections.
#
# Example:
# MIN_INTERSECTION_SIZE = 1000
MIN_INTERSECTION_SIZE = None


# Five sets have a maximum of 31 non-empty combinations.
MAX_INTERSECTIONS = 31


# Rotation angle for intersection count labels.
LABEL_ROTATION = 45


# Figure dimensions in inches.
FIGURE_WIDTH = 16
FIGURE_HEIGHT = 8


# ============================================================
# Helper functions
# ============================================================

def quote_identifier(identifier):
    """
    Safely quote a DuckDB table or column identifier.

    This is necessary for column names such as Arab-PG because
    DuckDB would otherwise interpret the hyphen as subtraction.
    """
    return '"' + identifier.replace('"', '""') + '"'


def format_count(value, _position=None):
    """
    Format large values using compact suffixes.

    Examples:
        2,604,236,942 -> 2.60B
        687,297,066   -> 687.3M
        97,159        -> 97.2K
        300           -> 300
    """
    value = float(value)

    if abs(value) >= 1_000_000_000:
        return "{:.2f}B".format(
            value / 1_000_000_000
        )

    if abs(value) >= 1_000_000:
        return "{:.1f}M".format(
            value / 1_000_000
        )

    if abs(value) >= 1_000:
        return "{:.1f}K".format(
            value / 1_000
        )

    return "{:,.0f}".format(value)


def check_table_exists(connection):
    """Confirm that the requested table exists."""
    result = connection.execute(
        """
        SELECT COUNT(*)
        FROM information_schema.tables
        WHERE table_name = ?
        """,
        [TABLE_NAME],
    ).fetchone()

    table_exists = result[0]

    if table_exists:
        return

    available_tables = connection.execute(
        """
        SELECT
            table_schema,
            table_name
        FROM information_schema.tables
        ORDER BY
            table_schema,
            table_name
        """
    ).fetchdf()

    raise RuntimeError(
        "Table '{}' was not found.\n\n"
        "Available tables:\n{}".format(
            TABLE_NAME,
            available_tables.to_string(index=False),
        )
    )


def check_required_columns(connection):
    """Confirm that all required columns are present."""
    quoted_table = quote_identifier(TABLE_NAME)

    describe_rows = connection.execute(
        "DESCRIBE {}".format(quoted_table)
    ).fetchall()

    available_columns = {
        row[0]
        for row in describe_rows
    }

    required_columns = set(
        ["kmer"] + list(SETS.values())
    )

    missing_columns = (
        required_columns - available_columns
    )

    if missing_columns:
        raise RuntimeError(
            "The following required columns are missing:\n"
            "  {}\n\n"
            "Available columns:\n"
            "  {}".format(
                sorted(missing_columns),
                sorted(available_columns),
            )
        )


def aggregate_intersections(connection):
    """
    Aggregate exact presence/absence combinations in DuckDB.

    For five binary sets, this returns at most 32 rows rather than
    transferring every k-mer into pandas.
    """
    select_expressions = []
    group_expressions = []

    for display_name, database_column in SETS.items():
        quoted_column = quote_identifier(
            database_column
        )

        quoted_alias = quote_identifier(
            display_name
        )

        # Any non-zero value is treated as present.
        # NULL values are treated as absent.
        expression = (
            "COALESCE({}, 0) <> 0".format(
                quoted_column
            )
        )

        select_expressions.append(
            "{} AS {}".format(
                expression,
                quoted_alias,
            )
        )

        group_expressions.append(expression)

    query = """
        SELECT
            {select_columns},
            COUNT(*)::UBIGINT AS intersection_size
        FROM {table_name}
        GROUP BY
            {group_columns}
        ORDER BY
            intersection_size DESC
    """.format(
        select_columns=", ".join(
            select_expressions
        ),
        table_name=quote_identifier(
            TABLE_NAME
        ),
        group_columns=", ".join(
            group_expressions
        ),
    )

    print(
        "Aggregating exact intersections "
        "inside DuckDB..."
    )

    return connection.execute(query).fetchdf()


def prepare_intersections(
    intersections,
):
    """
    Prepare aggregated intersection data for upsetplot.

    Returns
    -------
    Tuple[pandas.DataFrame, pandas.Series]
        The human-readable intersection table and the MultiIndex
        Series required by upsetplot.
    """
    if intersections.empty:
        raise RuntimeError(
            "No rows were returned from the "
            "k-mer presence table."
        )

    set_names = list(SETS.keys())

    for column in set_names:
        intersections[column] = (
            intersections[column].astype(bool)
        )

    intersections["intersection_size"] = (
        intersections[
            "intersection_size"
        ].astype("uint64")
    )

    # Remove rows absent from every reference.
    present_in_at_least_one = (
        intersections[set_names].any(axis=1)
    )

    all_absent_count = intersections.loc[
        ~present_in_at_least_one,
        "intersection_size",
    ].sum()

    if all_absent_count > 0:
        print(
            "Removing all-absent records: "
            "{:,}".format(
                int(all_absent_count)
            )
        )

    intersections = intersections.loc[
        present_in_at_least_one
    ].copy()

    if MIN_INTERSECTION_SIZE is not None:
        intersections = intersections.loc[
            intersections["intersection_size"]
            >= MIN_INTERSECTION_SIZE
        ].copy()

    if intersections.empty:
        raise RuntimeError(
            "No intersections remain after filtering. "
            "Reduce MIN_INTERSECTION_SIZE."
        )

    intersections["degree"] = (
        intersections[set_names].sum(axis=1)
    )

    intersections = intersections.sort_values(
        by=[
            "intersection_size",
            "degree",
        ],
        ascending=[
            False,
            False,
        ],
    )

    upset_series = intersections.set_index(
        set_names
    )["intersection_size"]

    return intersections, upset_series


def add_intersection_labels(axis):
    """
    Add compact, angled labels above the vertical bars.
    """
    for bar in axis.patches:
        height = bar.get_height()

        if height <= 0:
            continue

        x_position = (
            bar.get_x()
            + bar.get_width() / 2
        )

        axis.annotate(
            format_count(height),
            xy=(
                x_position,
                height,
            ),
            xytext=(2, 5),
            textcoords="offset points",
            ha="left",
            va="bottom",
            rotation=LABEL_ROTATION,
            rotation_mode="anchor",
            fontsize=7,
            clip_on=False,
        )


def add_total_labels(axis):
    """
    Add compact labels beside the horizontal total bars.
    """
    for bar in axis.patches:
        width = abs(bar.get_width())

        if width <= 0:
            continue

        left_edge = min(
            bar.get_x(),
            bar.get_x() + bar.get_width(),
        )

        y_position = (
            bar.get_y()
            + bar.get_height() / 2
        )

        axis.annotate(
            format_count(width),
            xy=(
                left_edge,
                y_position,
            ),
            xytext=(-5, 0),
            textcoords="offset points",
            ha="right",
            va="center",
            fontsize=7,
            clip_on=False,
        )


def create_upset_plot(upset_series):
    """Create and save the UpSet plot."""
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    with warnings.catch_warnings():
        warnings.simplefilter(
            "ignore",
            category=FutureWarning,
        )

        upset = UpSet(
            upset_series,
            subset_size="sum",
            sort_by="cardinality",
            sort_categories_by="cardinality",

            # Disable automatic count labels because they overlap.
            show_counts=False,

            show_percentages=False,
            element_size=38,
            intersection_plot_elements=8,
            totals_plot_elements=3,
            max_subset_rank=MAX_INTERSECTIONS,
        )

        axes = upset.plot()

    figure = plt.gcf()

    figure.set_size_inches(
        FIGURE_WIDTH,
        FIGURE_HEIGHT,
    )

    intersection_axis = axes[
        "intersections"
    ]

    totals_axis = axes[
        "totals"
    ]

    # Add manually formatted labels.
    add_intersection_labels(
        intersection_axis
    )

    add_total_labels(
        totals_axis
    )

    # Titles and axis labels.
    intersection_axis.set_title(
        "Exact k-mer intersections across "
        "pangenome references",
        pad=38,
    )

    intersection_axis.set_ylabel(
        "Number of k-mers"
    )

    totals_axis.set_xlabel(
        "Total k-mers per reference"
    )

    # Compact tick labels.
    intersection_axis.yaxis.set_major_formatter(
        plt.FuncFormatter(
            format_count
        )
    )

    totals_axis.xaxis.set_major_formatter(
        plt.FuncFormatter(
            format_count
        )
    )

    # Add space above bars for angled labels.
    lower_limit, upper_limit = (
        intersection_axis.get_ylim()
    )

    intersection_axis.set_ylim(
        lower_limit,
        upper_limit * 1.16,
    )

    # Remove unnecessary plot borders.
    for axis in (
        intersection_axis,
        totals_axis,
    ):
        axis.spines[
            "top"
        ].set_visible(False)

        axis.spines[
            "right"
        ].set_visible(False)

    figure.suptitle(
        "Pangenome k-mer overlap",
        fontsize=14,
        fontweight="bold",
        y=1.01,
    )

    figure.savefig(
        OUTPUT_PNG,
        dpi=400,
        bbox_inches="tight",
    )

    figure.savefig(
        OUTPUT_PDF,
        bbox_inches="tight",
    )

    plt.close(figure)


# ============================================================
# Main
# ============================================================

def main():
    """Run the complete k-mer UpSet analysis."""
    OUTPUT_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not DATABASE.exists():
        raise FileNotFoundError(
            "DuckDB database not found:\n"
            "{}".format(DATABASE)
        )

    print("=" * 70)
    print("K-mer UpSet plot")
    print("=" * 70)
    print(
        "Database: {}".format(
            DATABASE
        )
    )
    print(
        "Table:    {}".format(
            TABLE_NAME
        )
    )
    print()

    connection = duckdb.connect(
        str(DATABASE),
        read_only=True,
    )

    try:
        check_table_exists(
            connection
        )

        check_required_columns(
            connection
        )

        intersections = (
            aggregate_intersections(
                connection
            )
        )

    finally:
        connection.close()

    intersections, upset_series = (
        prepare_intersections(
            intersections
        )
    )

    intersections.to_csv(
        OUTPUT_COUNTS,
        index=False,
    )

    total_kmers = int(
        upset_series.sum()
    )

    print()
    print("Intersection counts:")
    print(
        intersections.to_string(
            index=False
        )
    )

    print()
    print(
        "Total represented k-mers: "
        "{:,}".format(total_kmers)
    )

    print(
        "Observed non-empty intersections: "
        "{:,}".format(
            len(upset_series)
        )
    )

    print()
    print(
        "Creating UpSet plot..."
    )

    create_upset_plot(
        upset_series
    )

    print()
    print("Outputs written:")

    print(
        "  PNG:    {}".format(
            OUTPUT_PNG
        )
    )

    print(
        "  PDF:    {}".format(
            OUTPUT_PDF
        )
    )

    print(
        "  Counts: {}".format(
            OUTPUT_COUNTS
        )
    )

    print("=" * 70)


if __name__ == "__main__":
    main()
