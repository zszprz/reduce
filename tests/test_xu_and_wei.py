import pytest
import numpy as np
from reduce import *

def test_xu_and_wei_cr():
    # Create a random 4x4 matrix
    matrix = np.random.rand(4, 4)

    # Calculate the initial Consistency Ratio (CR)
    matrix_size = matrix.shape[0]
    ci_initial = (return_max_eigenvalue(matrix) - matrix_size) / (matrix_size - 1)
    ri_initial = return_ri(matrix_size)
    cr_initial = ci_initial / ri_initial

    # Apply the xu_and_wei_cr function
    lambd = 0.5
    threshold = 0.1
    new_matrix, ci_final, cr_final = xu_and_wei_cr(matrix, lambd, threshold)

    # Assert that the new CR is less than the threshold
    assert cr_final < threshold

    # Assert that the new CR is less than the initial CR
    assert cr_final < cr_initial

    # Assert that the new matrix has the same size as the original matrix
    assert new_matrix.shape == matrix.shape
