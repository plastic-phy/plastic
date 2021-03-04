# distutils: sources = silly-test.c

cimport sillytestpyspecs as stp

def add_and_make_sabbia(fst, snd):

    trd = stp.add(fst, snd)
    cdef stp.sabbia* sabbia = stp.make_sabbia()

    return trd, sabbia.el

