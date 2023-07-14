import pytest
import numpy as np
from reduce import *

def test_create_pc_matrices():
    matrix = create_pc_matrices(3)
    assert matrix.shape == (3, 3)

def test_return_random_number():
    number = return_random_number()
    assert number >= 1/9 and number <= 9

def test_return_max_eigenvalue():
    matrix = np.array([[1, 2], [3, 4]])
    max_eigenvalue = return_max_eigenvalue(matrix)
    assert max_eigenvalue == pytest.approx(5.3722813232690143)

def test_calc_vecs():
    matrix = np.array([[1, 2], [3, 4]])
    vecs = calc_vecs(matrix)
    assert len(vecs) == 2

def test_calc_geo_mean():
    matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    geo_mean = calc_geo_mean(matrix)
    assert len(geo_mean) == 3

def test_return_ri():
    ri = return_ri(3)
    assert ri == 0.5247

def test_koczkodaj_index():
    matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    ki = koczkodaj_index(matrix)
    assert ki >= 0

def test_golden_wang_index():
    matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    gwi = golden_wang_index(matrix)
    assert gwi >= 0

def test_pelaez_lamata_index():
    matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    pli = pelaez_lamata_index(matrix)
    assert pli >= 0

def test_geometric_consistency_index():
    matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    gci = geometric_consistency_index(matrix)
    assert gci >= 0

def test_triads_geometric_consistency_index():
    matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    tgci = triads_geometric_consistency_index(matrix)
    assert tgci >= 0

def test_relative_error_index():
    matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    rei = relative_error_index(matrix)
    assert rei >= 0

def test_harmonic_consistency_index():
    matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    hci = harmonic_consistency_index(matrix)
    assert hci >= 0
