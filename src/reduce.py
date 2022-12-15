import random
import numpy as np
import scipy.linalg as la
from scipy.stats.mstats import gmean
from sympy import *
import csv
import math
from operator import itemgetter


# support functions

# function that create random pairwise comparison matrix

def create_pc_matrices(size):
    n_matrix = eye(size)  # generates identity matrix in size n x n

    for x in range(size):  # replace all zeros in that identity matrix with random picked numbers from defined list
        for y in range(size):
            if (n_matrix[x, y]) == 0:
                n_matrix[x, y] = return_random_number()

    for x in range(
            size):  # converts all elements in the lower half of matrix to reciprocals
        for y in range(x + 1):
            n_matrix[x, y] = 1 / n_matrix[y, x]

    return n_matrix


# function that return random number according to construction of pairwise comparison matrices

def return_random_number():
    choice_list = [Rational(1, 9), Rational(1, 8), Rational(1, 7), Rational(1, 6), Rational(1, 5), Rational(1, 4),
                   Rational(1, 3), Rational(1, 2), 1, 2, 3, 4, 5, 6, 7, 8, 9]

    return random.choice(choice_list)


# function that return max eigenvalue needed in further computations

def return_max_eigenvalue(n_matrix):
    nump_table = np.array(n_matrix).astype(np.float64)
    evals, evecs = la.eig(nump_table)

    evals = evals.real
    eigen_value_max = max(evals)

    return eigen_value_max


# function that calculates eigenvectors for matrix

def calc_vecs(n_matrix):
    nump_table = np.array(n_matrix).astype(np.float64)
    evals, evecs = la.eig(nump_table)
    abs_real_evecs = np.absolute(evecs.real)

    sum_of_vec = 0.0
    sum_of_all_vec = 0.0

    for evec in abs_real_evecs:  # dynamic calculation of the sum of vectors
        sum_of_vec = sum_of_vec + evec

    for sevec in sum_of_vec:  # calculation of the sum of all vectors
        sum_of_all_vec = sum_of_all_vec + sevec

    i = 0
    list_of_vecs = []

    while i < len(abs_real_evecs):  # dynamic calculation of the list of vectors for and vectors
        list_of_vecs.append(abs_real_evecs[i][0] / sum_of_vec[i])
        i += 1

    return list_of_vecs


# function that calculate geometric mean for given matrix

def calc_geo_mean(n_matrix):
    nump_table = np.array(n_matrix).astype(np.float64)
    list_of_geo_mean = []
    matrix_size = nump_table.shape[0]

    for n in range(matrix_size):
        w = gmean(nump_table[n, :])
        list_of_geo_mean.append(w)

    return list_of_geo_mean


# function that calculate eij for given matrix

def calc_eij(n_matrix):
    tmp_mac = np.array(n_matrix).astype(np.float64)
    l_vector = calc_geo_mean(tmp_mac)
    matrix_size = tmp_mac.shape[0]
    list_of_e = eye(matrix_size)

    for i in range(matrix_size):
        for j in range(matrix_size):

            if i < j:
                # calculation of the sum of eij
                list_of_e[i, j] = np.abs(np.log(tmp_mac[i, j] * (l_vector[j] / l_vector[i])))
            else:
                list_of_e[i, j] = 0

    # returns the entire matrix eij
    return list_of_e


# function that return precalculated RI value for given matrix size

def return_ri(size):
    if size == 3:
        return 0.5247
    elif size == 4:
        return 0.8815
    elif size == 5:
        return 1.1086
    elif size == 6:
        return 1.2479
    elif size == 7:
        return 1.3417
    elif size == 8:
        return 1.4056
    elif size == 9:
        return 1.4499
    elif size == 10:
        return 1.4854
    elif size == 11:
        return 1.5140
    elif size == 12:
        return 1.5365
    elif size == 13:
        return 1.5551
    elif size == 14:
        return 1.5713
    else:
        return 0


# function that import pairwise comparison matrices from prepared CSV file

def import_pc_matrix_from_csv(filename):
    reader = csv.reader(open(filename, "r"), delimiter=";")
    x = list(reader)
    imported_matrix = np.array(x).astype("float")

    matrix_size = len(x)
    matrix = eye(matrix_size)
    for n in range(matrix_size):
        for m in range(matrix_size):
            matrix[n, m] = imported_matrix[n][m]
    return matrix


# support function that convert given matrix to format needed for GMM index

def convert(imp_matrix):
    matrix = np.squeeze(np.asarray(imp_matrix))
    return matrix


# function that calculate dynamic vectors for given matrix

def dynamic_vectors(n_matrix):
    nump_table = np.array(n_matrix).astype(np.float64)
    evals, evecs = la.eig(nump_table)
    abs_real_evecs = np.absolute(evecs.real)
    vector_sum = 0
    list_of_vecs = []

    for i in range(len(n_matrix)):
        vector_sum = vector_sum + abs_real_evecs[i][0]

    for i in range(len(n_matrix)):
        list_of_vecs.append(abs_real_evecs[i][0] / vector_sum)

    return list_of_vecs


# function that calculate GMM vectors

def gmm_vectors(matrix):
    suma = 0
    geo_tmp = []
    tmp = 0
    vectors = []

    matrix_tmp = convert(matrix).astype('float64')
    size = int((len(matrix)))

    for i in range(size):
        geo = gmean(matrix_tmp[i], axis=0)
        geo_tmp.append(geo)

    for i in range(size):
        suma = suma + geo_tmp[i]

    for i in range(size):
        tmp = geo_tmp[i] / suma
        vectors.append(tmp)

    return vectors


# CR reduction algorithms

# function that implement Xu and Wei CR reduction algorithm and returns new matrix, CI and CR for given matrix,
# lambda parameter and CR threshold

def xu_and_wei_cr(matrix, lambd, threshold):
    tmp_mac = matrix
    matrix_size = int(math.sqrt(len(tmp_mac)))

    ci = (return_max_eigenvalue(tmp_mac) - matrix_size) / (matrix_size - 1)
    ri = return_ri(matrix_size)
    cr = ci / ri
    z = 0

    while cr >= threshold:

        z = z + 1
        l_vector = calc_vecs(tmp_mac)

        for n in range(matrix_size):
            for m in range(matrix_size):
                tmp_mac[n, m] = (pow(tmp_mac[n, m], lambd)) * (pow((l_vector[n] / l_vector[m]), (1 - lambd)))

        ci = (return_max_eigenvalue(tmp_mac) - matrix_size) / (matrix_size - 1)
        cr = ci / ri

    return tmp_mac, ci, cr


# function that implement Cao CR reduction algorithm and returns new matrix, CI and CR for given matrix, lambda parameter and CR threshold

def cao_cr(matrix, lambd, threshold):
    tmp_mac = matrix
    matrix_size = int(math.sqrt(len(tmp_mac)))

    ci = (return_max_eigenvalue(tmp_mac) - matrix_size) / (matrix_size - 1)
    ri = return_ri(matrix_size)
    cr = ci / ri
    z = 0

    while cr >= threshold:

        z = z + 1
        l_vector = calc_vecs(tmp_mac)

        for n in range(matrix_size):
            for m in range(matrix_size):
                tmp_mac[n, m] = (l_vector[n] / l_vector[m]) * (
                            lambd * tmp_mac[n, m] * (l_vector[m] / l_vector[n]) + (1 - lambd))

        ci = (return_max_eigenvalue(tmp_mac) - matrix_size) / (matrix_size - 1)
        cr = ci / ri

    return tmp_mac, ci, cr


# function that implement Szybowski CR reduction algorithm and returns new matrix, CI and CR for given matrix and CR threshold

def szybowski_cr(matrix, threshold):

    tmp_mac = matrix
    matrix_size = int(math.sqrt(len(tmp_mac)))

    ci = (return_max_eigenvalue(tmp_mac) - matrix_size) / (matrix_size - 1)
    ri = return_ri(matrix_size)
    cr = ci / ri
    z = 0

    while cr > threshold:

        z = z + 1
        list_of_eij = calc_eij(tmp_mac)
        l_vector = calc_geo_mean(tmp_mac)
        max_value = 0

        p = 0
        q = 0
        i = 0
        j = 0

        for i in range(matrix_size):
            for j in range(matrix_size):
                if i < j:
                    if max_value < list_of_eij[i, j]:
                        max_value = list_of_eij[i, j]
                        p = i
                        q = j

        tmp_mac[p, q] = l_vector[p] / l_vector[q]
        tmp_mac[q, p] = tmp_mac[p, q] ** -1

        ci = (return_max_eigenvalue(tmp_mac) - matrix_size) / (matrix_size - 1)
        cr = ci / ri

    return tmp_mac, ci, cr


# different indexes for Pairwise Comparison Matrices

# function that return Koczkodaj Index (KI) for given matrix

def koczkodaj_index(matrix):
    min_list = []
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            for k in range(len(matrix)):
                min_a = abs(1 - matrix[i][j] / (matrix[i][k] * matrix[k][j]))
                min_b = abs(1 - (matrix[i][k] * matrix[k][j]) / matrix[i][j])
                min_list.append(min(min_a, min_b))

    del matrix

    return max(min_list)


# function that returns Golden Wang Index (GWI) for given matrix

def golden_wang_index(matrix):
    sum_gw = []
    gw_index = 0

    for n in range(len(matrix)):
        sum_gw = 0
        for m in range(len(matrix)):
            sum_gw = sum_gw + matrix[m][n]
        sum_gw.append(sum_gw)

    for n in range(len(matrix)):
        for m in range(len(matrix)):
            matrix[m][n] = matrix[m][n] / sum_gw[n]

    vectors = dynamic_vectors(matrix)

    for n in range(len(matrix)):
        for m in range(len(matrix)):
            gw_index = gw_index + abs(matrix[n][m] - vectors[n])

    gw_index = (1 / (len(matrix) * len(matrix))) * gw_index

    del matrix
    return gw_index


# function that returns Palaez Lamata Index (PLI) for given matrix

def pelaez_lamata_index(matrix):
    pli = 0

    for i in range(len(matrix) - 2):
        for j in range(len(matrix) - 1):
            for k in range(len(matrix)):
                if j > i:
                    if k > j:
                        pli = pli + ((matrix[i][k] / (matrix[i][j] * matrix[j][k])) + (
                                    (matrix[i][j] * matrix[j][k]) / matrix[i][k]) - 2)

    n = len(matrix)
    pli = pli * (6 / (n * (n - 1) * (n - 2)))

    del matrix
    return pli


# function that returns Geometric Consistency Index (GCI) for given matrix

def geometric_consistency_index(matrix):
    gci = 0
    vectors = gmm_vectors(matrix)

    for i in range(len(matrix) - 1):
        for j in range(len(matrix)):
            if j > i:
                gci = gci + (log((matrix[i][j]) * (vectors[j] / vectors[i]))) * (
                    log((matrix[i][j]) * (vectors[j] / vectors[i])))

    n = len(matrix)
    gci = gci * (2 / ((n - 1) * (n - 2)))

    del matrix
    return gci


# function that returns Triads Geometric Consistency Index (TGCI) for given matrix

def triads_geometric_consistency_index(matrix):
    tgci = 0
    n = len(matrix)

    for k in range(len(matrix)):
        for j in range(len(matrix)):
            for i in range(len(matrix)):
                if k > j:
                    if j > i:
                        tgci = tgci + (log(matrix[i][j] * matrix[j][k] * matrix[k][i])) * (
                            log(matrix[i][j] * matrix[j][k] * matrix[k][i]))

    tgci = 2 * tgci
    tgci = tgci / (n * (n - 1) * (n - 2))

    del matrix
    return tgci


# function that returns Relative Error Index (REI) for given matrix

def relative_error_index(matrix):
    b = 0
    a = 0
    rei = 0
    vectors = gmm_vectors(matrix)

    for i in range(len(matrix)):
        for j in range(len(matrix)):
            a = a + (log(matrix[i][j]) - log((vectors[i] / vectors[j]))) * (
                        log(matrix[i][j]) - log((vectors[i] / vectors[j])))
            b = b + log(matrix[i][j]) * log(matrix[i][j])
    rei = a / b

    del a
    del b
    del matrix

    return rei


# function that returns Harmonic Consistency Index (HCI) for given matrix

def harmonic_consistency_index(matrix):
    sj = []

    for j in range(len(matrix)):
        tmp = 0
        for i in range(len(matrix)):
            tmp = matrix[i][j] + tmp
        sj.append(tmp)

    hm = 0

    for j in range(len(matrix)):
        hm = hm + (1 / sj[j])
    hm = (len(matrix)) / hm

    hci = ((hm - len(matrix)) * (len(matrix) + 1)) / (len(matrix) * (len(matrix) - 1))

    del matrix
    return hci
