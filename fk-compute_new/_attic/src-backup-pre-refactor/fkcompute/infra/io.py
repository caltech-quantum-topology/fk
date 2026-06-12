"""
File I/O utilities for FK computation.

This module provides utility functions for file operations
used throughout the FK computation package.
"""

import csv
import gzip
from typing import Any, Dict, List


def sort_any(xs: List[Any]) -> List[Any]:
    """
    Sort any list by string representation of elements.

    Parameters
    ----------
    xs
        List to sort.

    Returns
    -------
    list
        Sorted list.
    """
    return list(sorted(xs, key=lambda x: str(x)))


def find_where(lst: List[Any], predicate) -> Any:
    """
    Find the first element in a list satisfying a predicate.

    Parameters
    ----------
    lst
        List to search.
    predicate
        Function returning True for the desired element.

    Returns
    -------
    Any
        First matching element, or None if not found.
    """
    for x in lst:
        if predicate(x):
            return x
    return None


def find_index(lst: List[Any], predicate) -> int:
    """
    Find the index of the first element satisfying a predicate.

    Parameters
    ----------
    lst
        List to search.
    predicate
        Function returning True for the desired element.

    Returns
    -------
    int
        Index of first matching element, or None if not found.
    """
    for i, x in enumerate(lst):
        if predicate(x):
            return i
    return None


def csv_to_dicts(csv_file_path: str) -> List[Dict[str, str]]:
    """
    Read a CSV file into a list of dictionaries.

    Supports both plain CSV and gzip-compressed files.

    Parameters
    ----------
    csv_file_path
        Path to the CSV file (may end with .gz).

    Returns
    -------
    list[dict]
        List of dictionaries, one per row.
    """
    data_list = []

    # Open file as either gzipped or plain text
    open_func = gzip.open if csv_file_path.endswith('.gz') else open

    with open_func(csv_file_path, 'rt') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            data_list.append(row)

    return data_list


def tsv_to_dicts(tsv_file_path: str) -> List[Dict[str, str]]:
    """
    Read a TSV file into a list of dictionaries.

    Parameters
    ----------
    tsv_file_path
        Path to the TSV file.

    Returns
    -------
    list[dict]
        List of dictionaries, one per row.
    """
    data_list = []
    with open(tsv_file_path, 'r') as tsv_file:
        csv_reader = csv.DictReader(tsv_file, delimiter='\t')
        for row in csv_reader:
            data_list.append(row)
    return data_list


def save_dicts_to_tsv(data: List[Dict[str, Any]], file_path: str) -> None:
    """
    Save a list of dictionaries to a TSV file.

    Parameters
    ----------
    data
        List of dictionaries to save.
    file_path
        Path to the output TSV file.
    """
    with open(file_path, 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=data[0].keys(), delimiter='\t')
        writer.writeheader()
        for row in data:
            writer.writerow(row)
