from matrix_multiplication import multiply_single_thread, multiply_mp, multiply_mt, \
    multiply_cython, multiply_opencl, mse

EPS = 0.0000001

M_small = (2, 3)
N_small = (3, 2)

size = 100

M_huge = (size, size)
N_huge = (size, size)

RUNS = 10


def measure_time(func):
    import time
    import numpy as np
    runs = []

    for i in range(RUNS):
        start_time = time.time()
        result = func()
        runs.append((time.time() - start_time) * 1000)
    return np.mean(runs)


def run_test(M, N, multiply_func, request):
    import numpy as np

    a = np.random.rand(M[0], M[1]).astype(np.float32)
    b = np.random.rand(N[0], N[1]).astype(np.float32)

    numpy_result = np.dot(a, b)
    print()
    print("> %s" % request.node.name)
    print("Average (of %d runs): %.3f ms" % (RUNS, measure_time(lambda *_: multiply_func(a, b))))
    print()


# ===============

def test_time_small_matrix_single_thread(request):
    run_test(M_small, N_small, multiply_single_thread, request)


def test_time_huge_matrix_single_thread(request):
    run_test(M_huge, N_huge, multiply_single_thread, request)


# ===============

def test_time_small_matrix_multi_processing(request):
    run_test(M_small, N_small, multiply_mp, request)


def test_time_huge_matrix_multi_processing(request):
    run_test(M_huge, N_huge, multiply_mp, request)


# ===============

def test_time_small_matrix_multi_threading(request):
    run_test(M_small, N_small, multiply_mt, request)


def test_time_huge_matrix_multi_threading(request):
    run_test(M_huge, N_huge, multiply_mt, request)


# ===============

def test_time_small_matrix_cython(request):
    run_test(M_small, N_small, multiply_cython, request)


def test_time_huge_matrix_cython(request):
    run_test(M_huge, N_huge, multiply_cython, request)


# ===============

def test_time_small_matrix_opencl(request):
    run_test(M_small, N_small, multiply_opencl, request)


def test_time_huge_matrix_opencl(request):
    run_test(M_huge, N_huge, multiply_opencl, request)


# ===============

def check_correctness(M, N):
    import numpy as np
    a = np.random.rand(M[0], M[1]).astype(np.float32)
    b = np.random.rand(N[0], N[1]).astype(np.float32)

    np_result = np.dot(a, b)
    mp_result = multiply_mp(a, b)
    mt_result = multiply_mt(a, b)
    cython_result = multiply_cython(a, b)
    cl_result = multiply_opencl(a, b)

    assert mse(np_result, mp_result) < EPS
    assert mse(np_result, mt_result) < EPS
    assert mse(np_result, cython_result) < EPS
    assert mse(np_result, cl_result) < EPS


def test_correctness():
    get_range = lambda *_: range(2, 5)
    for i in get_range():
        for j in get_range():
            for k in get_range():
                check_correctness((i, k), (k, j))
