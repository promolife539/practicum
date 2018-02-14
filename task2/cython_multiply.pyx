import numpy as np
cimport numpy as np

def multiply(np.ndarray a, np.ndarray b):
    cdef np.ndarray result = np.zeros((a.shape[0], b.shape[1]))

    cdef int i, j, k

    for i in range(a.shape[0]):
        for j in range(b.shape[1]):
            for k in range(a.shape[1]):
                result[i, j] += a[i, k] * b[k, j]

    return result
