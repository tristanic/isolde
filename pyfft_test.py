import numpy

from time import time
import pyfftw

def test_fft(api, data):
    if api == 'pyfftw':
        start_time = time()
        b = pyfftw.interfaces.numpy_fft.fftn(data)
        print 'Forward fft took ' + str(time() - start_time) + ' seconds.'
        start_time = time()
        b = pyfftw.interfaces.numpy_fft.fftn(data)
        print 'Repeat fft took ' + str(time() - start_time) + ' seconds.'
        start_time = time()
        c = pyfftw.interfaces.numpy_fft.ifftn(b)
        print 'Inverse fft took ' + str(time() - start_time) + ' seconds.'
        start_time = time()
        c = pyfftw.interfaces.numpy_fft.ifftn(b)
        print 'Repeat inverse fft took ' + str(time() - start_time) + ' seconds.'
        error = numpy.abs(numpy.sum(numpy.abs(data) - numpy.abs(c)))/data.size
        print 'Average round-trip error = ' + str(error)
        return
    
    elif api == 'CUDA':
        from pyfft.cuda import Plan
        import pycuda.driver as cuda
        from pycuda.tools import make_default_context
        import pycuda.gpuarray as gpuarray
    
        start_time = time()
        cuda.init()
        context = make_default_context()
        stream=cuda.Stream()
        print('CUDA initialisation took {} seconds'.format(time()-start_time))
        start_time=time()
        plan = Plan(data.shape, stream=stream)
        print('Plan creation and compiling took {} seconds'.format(time()-start_time))
        start_time = time()
        
        gpu_data = gpuarray.to_gpu(data)
        plan.execute(gpu_data)
        result = gpu_data.get()
        print('Forward FFT took {} seconds'.format(time()-start_time))
        start_time=time()
        plan.execute(gpu_data, inverse=True)
        result = gpu_data.get()
        print('Inverse FFT took {} seconds'.format(time()-start_time))
        error = numpy.abs(numpy.sum(numpy.abs(data) - numpy.abs(result)))/data.size
        print 'Average round-trip error = ' + str(error)
        context.pop()
        
    elif api == 'OpenCL':
        from pyfft.cl import Plan
        import pyopencl as cl
        import pyopencl.array as cl_array
        
        start_time = time()
        context = cl.create_some_context(interactive=False)
        queue = cl.CommandQueue(context)
        print('OpenCL initialisation took {} seconds'.format(time()-start_time))
        
        start_time=time()
        plan = Plan(data.shape, queue=queue)
        print('Plan creation and compiling took {} seconds'.format(time()-start_time))
        start_time = time()
        gpu_data=cl_array.to_device(queue, data)
        plan.execute(gpu_data.data)
        result = gpu_data.get()
        
        print('Forward FFT took {} seconds'.format(time()-start_time))
        start_time=time()
        plan.execute(gpu_data.data, inverse=True)
        result = gpu_data.get()
        print('Inverse FFT took {} seconds'.format(time()-start_time))
        error = numpy.abs(numpy.sum(numpy.abs(data) - numpy.abs(result)))/data.size
        print 'Average round-trip error = ' + str(error)
        
nx, ny, nz = (256,256,1024)
#nx, ny, nz = numpy.random.randint(0, 512, 3)
print 'Array dimensions are ' + ','.join((str(nx), str(ny), str(nz)))
data = numpy.random.rand(nx, ny, nz).astype(numpy.complex64)
a = pyfftw.empty_aligned((nx, ny, nz), dtype='complex64', n=16)
a[:] = data


#~ print 'pyfftw'
#~ test_fft('pyfftw', a)
print 'OpenCL'
test_fft('OpenCL', data) 
print 'CUDA'   
test_fft('CUDA', data)
