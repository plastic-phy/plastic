from pandas import DataFrame
import numpy as np
import tatsu as ts

class DuplicateLabelsError(Exception): pass
class EmptyMatrixError(Exception): pass
class MatrixLabelSizeMismatch(Exception): pass

class LabeledMutationMatrix:

    
    def __init__(self, mutation_matrix, cell_labels = None, mutation_labels = None):

        def is_collection(obj):
            try:
                list(obj)
                return True
            except TypeError:
                return False
        
        if not is_collection(mutation_matrix):
            raise TypeError('The matrix must be a collection, but {0} is not.'.format(mutation_matrix))

        height = len(mutation_matrix)
        if height == 0:
            raise EmptyMatrixError()

        for row in list(mutation_matrix):
            if not is_collection(row):
                raise TypeError('Each row of the matrix must be a collection, but {0} was not.'.format(row))
            for element in row:
                if element not in {1, 2, 0}:
                    raise ValueError(element)
        width = len(mutation_matrix[0])
        if any([len(row) != width for row in mutation_matrix]):
            raise NotAMatrixError('Each row in the matrix must be of the same length')

        if cell_labels is None:
            cell_labels = [str(i) for i in range(1, height + 1)]
        if mutation_labels is None:
            mutation_labels = [str(i) for i in range(1, width + 1)]

        def _validate_labels(label_list, expected_length):
            if not is_collection(label_list):
                raise TypeError('mutation_labels and cell_labels must be collections, but {0} is not'.format(label_list))
            if any([isinstance(label_list, tp) for tp in {str, bytes}]):
                raise TypeError('strings like {0} are not accepted as label lists'.format(label_list))
            
            if len(label_list) != expected_length:
                raise MatrixLabelSizeMismatch('expected {0} labels, but found {1} in {2}.'.format(expected_length, len(label_list), label_list))
            
            for label in label_list:
                if not isinstance(label, str):
                    raise TypeError('each label must be a string, but {0} is not.'.format(label))
                if len(label) == 0 or len(bytes(label, 'ascii')) > 254:
                    raise ValueError('The label "{0}" is either an empty string, or a label with more than 254 characters.'.format(label))

            if len(set(label_list)) != len(label_list):
                raise DuplicateLabelsError()

        _validate_labels(cell_labels, height)
        _validate_labels(mutation_labels, width)
        
        self._data = DataFrame(mutation_matrix, index = cell_labels, columns = mutation_labels, dtype = int, copy = True)

    # The initial choice is to make the matrix immutable and to use external representations that are as general
    # as possible.
        
    def matrix(self):
        return self._data.to_numpy().tolist()
    
    @property
    def cell_labels(self):
        return list(self._data.index)
    
    @property
    def mutation_labels(self):
        return list(self._data.columns)

    def cells(self):
        return {lb : list(self._data.loc[lb, :]) for lb in self.cell_labels}

    def mutations(self):
        return {lb : list(self._data.loc[:, lb]) for lb in self.mutation_labels}

    # Representing the serializable form with strings lets us use existing formats for the representations.
    def to_serializable_dict(self):
        lines = ["".join([" " + str(element) for element in line]).strip() + '\n' for line in self.matrix()]
        mutation_matrix_as_string = "".join(lines)
        out = {
            'mutation_matrix' : mutation_matrix_as_string,
            'cell_labels' : "".join([lb + '\n' for lb in self.cell_labels]),
            'mutation_labels' : "".join([lb + '\n' for lb in self.mutation_labels])
        }
        return out

    def _from_strings(mutation_matrix, cell_labels = None, mutation_labels = None, matstring_format = "SASC", **kwargs):
        
        def _parse_labels(labels_string):
            if labels_string is None:
                return None
            return [lb.strip() for lb in labels_string.splitlines() if lb.strip() != ""]
        
        return LabeledMutationMatrix(
            _MatrixParser(matstring_format).parse_matrix(mutation_matrix),
            _parse_labels(cell_labels),
            _parse_labels(mutation_labels)
        )

    def from_serializable_dict(dict_representation):
        return LabeledMutationMatrix._from_strings(matstring_format = 'SASC', **dict_representation)

    def from_files(matrix_file, matstring_format = 'SASC', cells_file = None, mutations_file = None):

        def _read_nullable(file_name, default_output = None):
            if file_name is None:
                return default_output
            else:
                with open(file_name, 'r') as f:
                    return f.read()

        with open(matrix_file, 'r') as f:
            mutation_matrix = f.read()

        return LabeledMutationMatrix._from_strings(
            mutation_matrix = mutation_matrix,
            cell_labels = _read_nullable(cells_file),
            mutation_labels = _read_nullable(mutations_file),
            matstring_format = matstring_format
        )

        
class _MatrixFileFormat:
    def __init__(self, grammar, value_map, transpose):
        self.grammar = grammar
        self.value_map = value_map
        self.transpose = transpose

# Should this go into its own file?
_formats = {
    'SASC' : _MatrixFileFormat(
        grammar = (
            r"""@@grammar :: SASC
            @@whitespace ::  /(?s)[ \t\r\f\v]+/
            start = matrix $ ;
            cell = "0" ~ | "1" ~ | "2" ~ ;
            row = {cell};
            matrix = ("\n").{ row } ;
            """
        ),
        value_map = {'0' : 0, '1': 1, '2' : 2},
        transpose = False
    ),
    
    'SCITE' : _MatrixFileFormat(
        grammar = (
            r"""
            @@grammar :: SCITE
            @@whitespace ::  /(?s)[ \t\r\f\v]+/
            start = matrix $ ;
            cell = "0" ~ | "1" ~ | "2" ~ | "3" ~ ;
            row = {cell};
            matrix = ("\n").{ row } ;
            """
        ),
        value_map = {'0' : 0, '1' : 1, '2' : 2, '3' : 2},
        transpose = False
    ),
    
    'SPHYR' : _MatrixFileFormat(
        grammar = (
            r"""
            @@grammar :: SPHYR
            @@whitespace ::  /(?s)[ \t\r\f\v]+/
            @@eol_comments :: /#([^\n]*?)$/
            start = cell_number_line "\n" snv_number_line "\n" matrix $ ;
            cell_number_line = /\d+/;
            snv_number_line = /\d+/;
            cell = "0" ~ | "1" ~ | "-1" ~ ;
            row = {cell};
            matrix = ("\n").{row};
            """
        ),
        value_map = {'0' : 0, '1': 1, '-1' : 2},
        transpose = True
    )
}


class NotAMatrixError(Exception): pass
class MeaninglessCellValueError(Exception): pass

class _MatrixParser:

    def __init__(self, parser_name):
        self.matrix_format = _formats[parser_name]
        self.inner_parser = ts.compile(self.matrix_format.grammar)

    def parse_matrix(self, matrix_string):
        parsed_string = [line for line in self.inner_parser.parse(matrix_string) if len(line) != 0]

        matrix = np.vectorize(self.matrix_format.value_map.get)(parsed_string)
        if self.matrix_format.transpose:
            matrix = np.transpose(matrix)
        return matrix

        
        
                

    

        
        

        



        