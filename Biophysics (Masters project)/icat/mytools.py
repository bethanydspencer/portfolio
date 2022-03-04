#!/usr/bin/env python
import numpy, array
from numpy import asarray
def read_array(filename, dtype, separator=' '):
    """ Read a file with an arbitrary number of columns.
        The type of data in each column is arbitrary
        It will be cast to the given dtype at runtime
    """
    cast = numpy.cast
    data = [[] for dummy in xrange(len(dtype))]
    for line in open(filename, 'r'):
        fields = line.strip().split(separator)
        for i, number in enumerate(fields):
            data[i].append(number)
    for i in xrange(len(dtype)):
        data[i] = cast[dtype[i]](data[i])
    return numpy.rec.array(data, dtype=dtype)


def myread(file, commentchar='#'):
    """Load a table with numbers into a two-dim. Numpy array."""
    # read until next blank line:
    r = []  # total set of numbers (r[i]: numbers in i-th row)
    while 1:    # might call read several times from a file
        line = file.readline()
        if not line: break  # end of file
        if line.isspace(): break# blank line
        if line[0] == commentchar: continue # treat next line
        r.append([float(s) for s in line.split()])
    return asarray(r, dtype=float)


def readfile(filename, commentchar='#'):
    """As read. Return columns as separate arrays."""
    f = open(filename, 'r')
    a = myread(f, commentchar)
    r = [a[:,i] for i in range(a.shape[1])]
    return r

def movie(outputFileName='output.avi', vcodec='mpeg4'):
    """
    String a bunch of png files together and create a movie.
    """
    #command to encode files into a movie
    command = 'mencoder mf://*.png -mf type=png:w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=' + vcodec + ' -oac copy -o ' + outputFileName
    
    from subprocess import call
    #Execute command
    retcode = call(command)