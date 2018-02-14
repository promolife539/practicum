import numpy as np
import multiprocessing as mp
import multiprocessing.dummy as mt
import cython_multiply
import opencl_multiply

def mse(a, b):
    return np.square(np.subtract(a, b)).mean()


def multiply_single_thread(a, b):
    result = np.zeros((a.shape[0], b.shape[1]))

    for i in range(a.shape[0]):
        for j in range(b.shape[1]):
            for k in range(a.shape[1]):
                result[i, j] += a[i, k] * b[k, j]

    return result


def mp_worker(args):
    i, a, b = args
    result = np.zeros(b.shape[1])
    for j in range(b.shape[1]):
        for k in range(a.shape[1]):
            result[j] += a[i, k] * b[k, j]
    return result


def multuply_multithreading(a, b, mp_lib):
    pool = mp_lib.Pool(mp.cpu_count())
    return np.array(pool.map(mp_worker, [(i, a, b) for i in range(a.shape[0])]))


def multiply_mp(a, b):
    return multuply_multithreading(a, b, mp)


def multiply_mt(a, b):
    return multuply_multithreading(a, b, mt)


def multiply_cython(a, b):
    return cython_multiply.multiply(a, b)


def multiply_opencl(a, b):
    return opencl_multiply.multiply(a, b)

### Testing

# a = np.random.rand(4, 7).astype(np.float32)
# b = np.random.rand(7, 4).astype(np.float32)
#
# np_result = np.dot(a, b)
# cl_result = multiply_opencl(a, b)
#
# print(np_result)
# print(cl_result)
# print(mse(np_result, cl_result))
