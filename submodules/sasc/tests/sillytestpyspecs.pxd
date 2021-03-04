cdef extern from "silly-test.h":

    struct s:
        int el
        int *arr

    ctypedef s sabbia

    sabbia *make_sabbia()
    int add(int fst, int snd)
