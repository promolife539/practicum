import numpy as np
import pyopencl as cl

# (MxN) * (NxK)
def multiply(a, b):
    m = a.shape[0]
    n = a.shape[1]
    k = b.shape[1]

    a_np = a.flatten()
    b_np = b.flatten()

    res_np = np.zeros(m * k).astype(np.float32)

    platforms = cl.get_platforms()
    ctx = cl.Context(
        dev_type=cl.device_type.ALL,
        properties=[(cl.context_properties.PLATFORM, platforms[0])])

    prg = cl.Program(ctx, """
        float get_val(__global const float *arr, int width, int y, int x);
        void add_val(__global float *arr, int width, int y, int x, float val);
        
        float get_val(__global const float *arr, int width, int y, int x) {
            return arr[y * width + x];
        }
        
        void add_val(__global float *arr, int width, int y, int x, float val) {
            arr[y * width + x] += val;
        }
        
        __kernel void multiply(
            __global const float *a_g, __global const float *b_g, __global float *res_g,
            const int m, const int n, const int k
        ){
            int i = get_global_id(0);

            for (int j = 0; j < k; j++) {
                for (int l = 0; l < n; l++) {
                    float val = get_val(a_g, n, i, l) * get_val(b_g, k, l, j);
                    add_val(res_g, k, i, j, val);
                } 
            }
        }
    """).build()

    mf = cl.mem_flags
    a_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=a_np)
    b_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=b_np)
    res_g = cl.Buffer(ctx, mf.WRITE_ONLY, res_np.nbytes)

    queue = cl.CommandQueue(ctx)
    prg.multiply(queue, [m], None, a_g, b_g, res_g, np.int32(m), np.int32(n), np.int32(k))
    cl.enqueue_copy(queue, res_np, res_g)

    return res_np.reshape((m, k))
