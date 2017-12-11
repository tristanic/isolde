import numpy
import reikna.cluda as cluda
from reikna.fft import FFT
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
        
    start_time = time()
    thr = api.Thread.create()
    print('Creating thread took {} seconds'.format(time() - start_time))
    start_time = time()
    fft = FFT(data, axes=(0,1,2))
    print('Preparing FFT took {} seconds'.format(time() - start_time))
    start_time = time()
    fftc = fft.compile(thr)
    print 'Compiling fft took ' + str(time() - start_time) + ' seconds.'
    
    start_time = time()
    data_dev = thr.to_device(data)
    fftc(data_dev, data_dev)
    transformed = data_dev.get()
    print 'Forward fft took ' + str (time() - start_time) + ' seconds.'
    
    start_time = time()
    fftc(data_dev, data_dev, inverse=True)
    back_transformed = data_dev.get()
    print 'Inverse fft took ' + str (time() - start_time) + ' seconds.'
    
    error = numpy.abs(numpy.sum(numpy.abs(data) - numpy.abs(back_transformed)))/data.size
    print 'Average round-trip error = ' + str(error)

nx= ny= nz = 509
#nx, ny, nz = numpy.random.randint(0, 512, 3)
print 'Array dimensions are ' + ','.join((str(nx), str(ny), str(nz)))
clapi = cluda.ocl_api()
cuapi = cluda.cuda_api()
data = numpy.random.rand(nx, ny, nz).astype(numpy.complex64)
a = pyfftw.empty_aligned((nx, ny, nz), dtype='complex64', n=16)
a[:] = data


print 'pyfftw'
test_fft('pyfftw', a)
print 'OpenCL'
test_fft(clapi, data) 
print 'CUDA'   
test_fft(cuapi, data)
