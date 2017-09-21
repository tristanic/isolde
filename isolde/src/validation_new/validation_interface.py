'''
All structural validation tasks should be run in (a) thread(s) to reduce their
impact on graphics performance. The interface between ISOLDE and the
validation thread(s) is defined here.
'''
import numpy
import multiprocessing as mp
from multiprocessing.pool import ThreadPool as Pool
