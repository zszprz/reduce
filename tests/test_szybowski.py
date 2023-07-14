import pytest
import numpy as np
from reduce import *

def test_szybowski_cr():
    # Create a random 4x4 matrix
    matrix = np.random.rand(4, 4)

    # Calculate the initial Consistency Ratio (CR)
    matrix_size = matrix.shape[0]
    ci_initial = (return_max_eigenvalue(matrix) - matrix_size) / (matrix_size - 1)
    ri_initial = return_ri(matrix_size)
    cr_initial = ci_initial / ri_initial

    # Apply the szybowski_cr function
    threshold = 0.1
    new_matrix, ci_final, cr_final = szybowski_cr(matrix, threshold)

    # Assert that the new CR is less than the threshold
    assert cr_final < threshold

    # Assert that the new CR is less than the initial CR
    assert cr_final < cr_initial

    # Assert that the new matrix has the same size as the original matrix
    assert new_matrix.shape == matrix.shape
