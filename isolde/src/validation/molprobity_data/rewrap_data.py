# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 11-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll



import os
import numpy
import glob

def rewrap_data(fname):
    '''
    Rewraps the angles in a MolProbity rotamer dataset from (0..360) to
    (-180..180) for consistency with the Ramachandran set.
    '''
    with open(fname, 'rt') as fin, open(os.path.join('rota_rewrapped', fname), 'wt') as fout:
        fout.write(fin.readline())
        dim_line = fin.readline()
        fout.write(dim_line)
        ndim = int(dim_line.split()[-1])
        fout.write(fin.readline())
        for i in range(ndim):
            line = fin.readline().split()
            if line[3] == '360.0':
                line[2] = '{:.1f}'.format((float(line[2])-180.0))
                line[3] = '{:.1f}'.format((float(line[3])-180.0))
            fout.write(' '.join(line)+'\n')
        fout.write(fin.readline())
        data = numpy.loadtxt(fin)
        data[:,0:ndim][data[:,0:ndim]>180] -=360
        idx = numpy.lexsort([data[:,i] for i in reversed(range(ndim))])
        data = data[idx]
        print('loaded data')
        for line in data:
            out_str = ' '.join(['{:.1f}'.format(float(d)) for d in line[0:ndim]]) + ' ' + '{:.15f}'.format(float(line[-1]))+'\n'
            fout.write(out_str)
            fout.flush()
        #~ for line in fin:
            #~ line_data = numpy.array([float(d) for d in line.split()])
            #~ line_data[line_data>180.0] -= 360.0
            #~ out_str = ' '.join(['{:.1f}'.format(d) for d in line_data[0:ndim]]) + ' ' + '{:.15f}'.format(line_data[-1])+'\n'
            #~ fout.write(out_str)

for f in glob.glob('rota*.data'):
    rewrap_data(f)
