# type: ignore[attr-defined]
"""Reduce is a Python package that helps with pairwise comparison matrices computations.

This package provides a variety of functions for creating and manipulating pairwise comparison matrices,
as well as computing various indices and performing CR reductions.

Modules:
    reduce: Contains the main functionality for creating and manipulating pairwise comparison matrices.

Functions in the reduce module:
    create_pc_matrices(size): Creates a random pairwise comparison matrix of a given size.
    return_random_number(): Returns a random number according to the construction of pairwise comparison matrices.
    return_max_eigenvalue(n_matrix): Returns the maximum eigenvalue needed in further computations.
    calc_vecs(n_matrix): Calculates eigenvectors for a given matrix.
    calc_geo_mean(n_matrix): Calculates the geometric mean for a given matrix.
    calc_eij(n_matrix): Calculates eij for a given matrix.
    return_ri(size): Returns the precalculated RI value for a given matrix size.
    import_pc_matrix_from_csv(filename): Imports a pairwise comparison matrix from a prepared CSV file.
    convert(imp_matrix): Converts a given matrix to the format needed for the GMM index.
    dynamic_vectors(n_matrix): Calculates dynamic vectors for a given matrix.
    gmm_vectors(matrix): Calculates GMM vectors.
    xu_and_wei_cr(matrix, lambd, threshold): Implements the Xu and Wei CR reduction algorithm.
    cao_cr(matrix, lambd, threshold): Implements the Cao CR reduction algorithm.
    szybowski_cr(matrix, threshold): Implements the Szybowski CR reduction algorithm.
    koczkodaj_index(matrix): Returns the Koczkodaj Index (KI) for a given matrix.
    golden_wang_index(matrix): Returns the Golden Wang Index (GWI) for a given matrix.
    pelaez_lamata_index(matrix): Returns the Palaez Lamata Index (PLI) for a given matrix.
    geometric_consistency_index(matrix): Returns the Geometric Consistency Index (GCI) for a given matrix.
    triads_geometric_consistency_index(matrix): Returns the Triads Geometric Consistency Index (TGCI) for a given matrix.
    relative_error_index(matrix): Returns the Relative Error Index (REI) for a given matrix.
    harmonic_consistency_index(matrix): Returns the Harmonic Consistency Index (HCI) for a given matrix.
"""

import sys

if sys.version_info >= (3, 8):
    from importlib import metadata as importlib_metadata
else:
    import importlib_metadata

from . import reduce

def get_version() -> str:
    """Return the current version of the reduce package."""
    try:
        return importlib_metadata.version(__name__)
    except importlib_metadata.PackageNotFoundError:  # pragma: no cover
        return "unknown"

version: str = get_version()
