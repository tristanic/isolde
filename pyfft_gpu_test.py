import pyfftw
import numpy
from time import time

a = pyfftw.empty_aligned((634,634,634), dtype='complex64', n=16)
a[:] = numpy.random.randn(634,634,634) + 1j*numpy.zeros((634,634,634))
start_time = time()
b = pyfftw.interfaces.numpy_fft.fft(a)
elapsed_time = time() - start_time
print 'FFTW 634^3 fft took ' + str(elapsed_time) + ' seconds.'

from pyfft.cuda import Plan
import numpy
import pycuda.driver as cuda
from pycuda.tools import make_default_context
import pycuda.gpuarray as gpuarray
from time import time


cuda.init()
context = make_default_context()
stream = cuda.Stream()

start_time = time()
plan = Plan((634, 634, 634), stream = stream)
elapsed_time = time() - start_time

print 'Generating CUDA plan took ' + str(elapsed_time) + ' seconds.'

data = (numpy.random.rand(634,634,634).astype(numpy.complex64) - 0.5 * 1000)
gpu_data = gpuarray.to_gpu(data)

start_time = time()
plan.execute(gpu_data)
result = gpu_data.get()
elapsed_time = time() - start_time

print '634^3 FFT on CUDA took ' + str(elapsed_time) + ' seconds.'

start_time = time()
plan.execute(gpu_data, inverse = True)
result = gpu_data.get()
elapsed_time = time() - start_time


print '634^3 inverse FFT on CUDA took ' + str(elapsed_time) + ' seconds.'

context.pop()


#####
# OpenCL
#####
from pyfft.cl import Plan
import pyopencl as cl
import pyopencl.array as cl_array

ctx = cl.create_some_context(interactive=False)
queue = cl.CommandQueue(ctx)

start_time = time()
plan=Plan((634,634,634), queue=queue)
elapsed_time = time() - start_time

print 'Generating OpenCL plan took ' + str(elapsed_time) + ' seconds.'


gpu_data = cl_array.to_device(queue, data)

start_time = time()
plan.execute(gpu_data.data)
result = gpu_data.get()
elapsed_time = time() - start_time

print '634^3 FFT on OpenCL took ' + str(elapsed_time) + ' seconds.'

start_time = time()
plan.execute(gpu_data.data, inverse = True)
result = gpu_data.get()
elapsed_time = time() - start_time


print '634^3 inverse FFT on OpenCL took ' + str(elapsed_time) + ' seconds.'




 
