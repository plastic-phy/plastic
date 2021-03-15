from labeledmutationmatrix import *
import pytest as pt
import pandas as pd
import numpy as np

def stub_matrix():
    return LabeledMutationMatrix([[0, 0, 0], [0, 0, 0], [1, 0, 0]], ['a', 'b', 'c'], ('sabbia', 'pollo', 'facocero'))

# TESTS ON INIT

def test_init():

    a = stub_matrix()
    assert a.matrix() == [[0, 0, 0], [0, 0, 0], [1, 0, 0]]
    assert a.cell_labels == ['a', 'b', 'c']
    assert a.mutation_labels == ['sabbia', 'pollo', 'facocero']

def test_init_without_explicit_labels():

    a = LabeledMutationMatrix([[0, 0, 0], [0, 0, 0]])

    assert a.matrix() == [[0, 0, 0], [0, 0, 0]]
    assert a.cell_labels == ['1', '2']
    assert a.mutation_labels == ['1', '2', '3']

def test_init_with_numpy_matrix():
    
    a = LabeledMutationMatrix(np.array([[0, 0, 0]]))

    assert a.matrix() == [[0, 0, 0]]

def test_init_with_doubles_in_matrix():
    a = LabeledMutationMatrix(np.array([[0.0, 0, 0]]))

    assert str(a.matrix()) == '[[0, 0, 0]]'

def test_notiterable_matrix_arg():

    with pt.raises(TypeError):
        LabeledMutationMatrix(None)

def test_not_two_dimensional_matrix_arg():

    with pt.raises(TypeError):
        LabeledMutationMatrix([2])

def test_ragged_matrix_arg():

    with pt.raises(NotAMatrixError):
        LabeledMutationMatrix([[0, 1, 2], [0, 1]])

def test_bad_iterable_arg():
    
    with pt.raises(ValueError):
        LabeledMutationMatrix("pollo")

def test_bad_values_in_matrix_arg():

    with pt.raises(ValueError):
        LabeledMutationMatrix([[69]])

def test_noniterable_label_arg():

    with pt.raises(TypeError):
        LabeledMutationMatrix([[0]], 3, ['a'])

def test_string_label_arg():

    with pt.raises(TypeError):
        LabeledMutationMatrix([[0]], ['1'], 'polloplex')

def test_wrong_length_label_arg():

    with pt.raises(MatrixLabelSizeMismatch):
        LabeledMutationMatrix([[0]], ['1', '2'], ['a'])

def test_duplicate_labels():

    with pt.raises(DuplicateLabelsError):
        LabeledMutationMatrix([[0, 0]], ['1'], ['1', '1'])

def test_bad_label():

    with pt.raises(TypeError):
        LabeledMutationMatrix([[0, 0]], [3], ['1', '1'])

def test_bad_length_label():

    with pt.raises(ValueError):
        LabeledMutationMatrix([[0]], ['a'*255], ['a'])

# CONVERSION_TO_SERALIZABLE_TESTS

def test_conversion_to_serializable():
    a = stub_matrix()
    d = a.to_serializable_dict()

    assert d['mutation_matrix'] == '0 0 0\n0 0 0\n1 0 0\n'
    assert d['cell_labels'] == 'a\nb\nc\n'
    assert d['mutation_labels'] == 'sabbia\npollo\nfacocero\n'

def test_conversion_roundtrip():
    a = stub_matrix()
    d = a.to_serializable_dict()
    a_after_rt = LabeledMutationMatrix.from_serializable_dict(d)

    assert a.matrix() == a_after_rt.matrix()
    assert a.cell_labels == a_after_rt.cell_labels
    assert a.mutation_labels == a_after_rt.mutation_labels

# CELLS AND MUTATIONS TESTS

def test_cells():

    a = stub_matrix()
    assert a.cells() == {'a' : [0, 0, 0], 'b' : [0, 0, 0], 'c' : [1, 0, 0]}
    assert isinstance(a.cells()['a'][0], int)

def test_mutations():

    a = stub_matrix()
    assert a.mutations() == {'sabbia': [0, 0, 1], 'pollo' : [0, 0, 0], 'facocero': [0, 0, 0]}

# IMMUTABILITY TESTS

def test_matrix_immutability():
    a = stub_matrix()
    a.matrix()[0][0] = 3
    b = stub_matrix()
    assert a.matrix() == b.matrix()

def test_labels_immutability():
    a = stub_matrix()
    a.cell_labels[0] = 'pollo'
    b = stub_matrix()

    assert a.cell_labels == b.cell_labels

# mutations uses the same logic, no need to test that as well
def test_cells_immutability():
    a = stub_matrix()
    a.cells()['a'][0] = 69
    b = stub_matrix()

    assert a.cells() == b.cells()





