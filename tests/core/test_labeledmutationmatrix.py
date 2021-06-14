from phylo._core.genotypematrix import *
import pytest as pt
import pandas as pd
import numpy as np

def stub_matrix():
    return GenotypeMatrix([[0, 0, 0], [0, 0, 0], [1, 0, 0]], ['a', 'b', 'c'], ('sabbia', 'pollo', 'facocero'))

class Test_Accessors:
    def test_cells(self):

        a = stub_matrix()
        assert a.cells() == {'a' : [0, 0, 0], 'b' : [0, 0, 0], 'c' : [1, 0, 0]}
        assert isinstance(a.cells()['a'][0], int)

    def test_mutations(self):

        a = stub_matrix()
        assert a.mutations() == {'sabbia': [0, 0, 1], 'pollo' : [0, 0, 0], 'facocero': [0, 0, 0]}


class Test_Immutability:

    def test_matrix_immutability_after_init_parameter_modification(self):
        mutable_init_matrix = np.array([[0, 0, 0], [1, 0, 1]])
        a = GenotypeMatrix(mutable_init_matrix)
        mutable_init_matrix[0, 0] = 69
        assert a.matrix()[0][0] == 0

    def test_labels_immutability_after_init_parameter_modification(self):
        mutable_init_labels = ['sabbia', 'polloplex']
        a = GenotypeMatrix([[0, 0]], mutation_labels = mutable_init_labels)
        mutable_init_labels[0] = 'you_thought_it_would_be_sabbia_but_it_is_me_dio'
        assert set(a.mutation_labels) == {'sabbia', 'polloplex'}

    def test_matrix_immutability_after_modification_to_accessor(self):
        a = stub_matrix()
        a.matrix()[0][0] = 3
        b = stub_matrix()
        assert np.all(a.matrix() == b.matrix())

    def test_labels_immutability_after_modification_through_accessor(self):
        a = stub_matrix()
        a.cell_labels[0] = 'pollo'
        b = stub_matrix()
        assert a.cell_labels == b.cell_labels

        # mutations uses the same logic, no need to test that as well
    def test_cells_immutability_after_modification_through_accessor(self):
        a = stub_matrix()
        a.cells()['a'][0] = 69
        b = stub_matrix()
        assert a.cells() == b.cells()


class TestInit:

    def test_init(self):

        a = stub_matrix()
        assert np.all(a.matrix() == np.array([[0, 0, 0], [0, 0, 0], [1, 0, 0]]))
        assert a.cell_labels == ['a', 'b', 'c']
        assert a.mutation_labels == ['sabbia', 'pollo', 'facocero']

    def test_init_without_explicit_labels(self):

        a = GenotypeMatrix([[0, 0, 0], [0, 0, 0]])
        assert np.all(a.matrix() == np.array([[0, 0, 0], [0, 0, 0]]))
        assert a.cell_labels == ['1', '2']
        assert a.mutation_labels == ['1', '2', '3']

    def test_init_with_numpy_matrix(self):

        a = GenotypeMatrix(np.array([[0, 0, 0]]))
        assert np.all(a.matrix() == np.array([[0, 0, 0]]))

    def test_init_with_doubles_in_matrix(self):

        a =GenotypeMatrix(np.array([[0.0, 0, 0]]))
        assert str(a.matrix()) == '[[0 0 0]]'

    def test_noncollection_matrix_arg(self):

        with pt.raises(TypeError):
            GenotypeMatrix(None)

    def test_not_two_dimensional_matrix_arg(self):

        with pt.raises(ValueError):
            GenotypeMatrix([2])

    def test_ragged_matrix_arg(self):

        with pt.raises(ValueError):
            GenotypeMatrix([[0, 1, 2], [0, 1]])

    def test_bad_collection_arg(self):

        with pt.raises(ValueError):
            GenotypeMatrix("pollo")

    def test_bad_values_in_matrix_arg(self):

        with pt.raises(ValueError):
            GenotypeMatrix([[69]])

    def test_noncollection_label_arg(self):

        with pt.raises(TypeError):
            GenotypeMatrix([[0]], 3, ['a'])

    def test_string_label_arg(self):

        with pt.raises(TypeError):
            GenotypeMatrix([[0]], ['1'], 'polloplex')

    def test_wrong_length_label_arg(self):

        with pt.raises(MatrixLabelSizeMismatch):
            GenotypeMatrix([[0]], ['1', '2'], ['a'])

    def test_duplicate_labels(self):

        with pt.raises(DuplicateLabelsError):
            GenotypeMatrix([[0, 0]], ['1'], ['1', '1'])

    def test_bad_label(self):

        with pt.raises(TypeError):
            GenotypeMatrix([[0, 0]], [3], ['1', '1'])

    def test_label_with_bad_characters(self):

        with pt.raises(ValueError):
            GenotypeMatrix([[0]], ["the issue with Ch- \nhey, we don't do this here!"], ['a'])


class Test_Serialization:
    def test_conversion_to_serializable(self):
        a = stub_matrix()
        d = a.to_serializable_dict()

        assert d['genotype_matrix'] == '0 0 0\n0 0 0\n1 0 0'
        assert d['cell_labels'] == 'a\nb\nc'
        assert d['mutation_labels'] == 'sabbia\npollo\nfacocero'

    def test_conversion_roundtrip(self):
        a = stub_matrix()
        d = a.to_serializable_dict()
        a_after_rt = GenotypeMatrix.from_serializable_dict(d)

        assert np.all(a.matrix() == a_after_rt.matrix())
        assert a.cell_labels == a_after_rt.cell_labels
        assert a.mutation_labels == a_after_rt.mutation_labels


class TestMatrixLoading:

    def test_plain_matrix_loading(self):

        sasc_format_matrix = GenotypeMatrix.from_files('matrepr_files/sasc_example')
