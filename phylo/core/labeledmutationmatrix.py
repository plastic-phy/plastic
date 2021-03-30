from pandas import DataFrame
import numpy as np
import tatsu as ts
import re

class DuplicateLabelsError(Exception): pass
class EmptyMatrixError(Exception): pass
class MatrixLabelSizeMismatch(Exception): pass

class LabeledMutationMatrix:

    
    def __init__(self, mutation_matrix, cell_labels = None, mutation_labels = None):
        """
        Validates a matrix that describes the results of an analysis of the mutations that have
        occurred in a sample of cells.

        Parameters:
        - mutation_matrix: A twice-subscriptable collection of values that can be converted to a
          bidimensional, nonempty np.ndarray of integer in which each value is either 0, 1 or 2.
        - cell_labels (list<string>?), by default None: Either None, or a list of strings with the
          same size as the height of the mutation_matrix representation as a np.ndarray. If the list exists,
          each element must be a non-empty string without spaces, newlines or other tabulation characters.
        - mutation_labels (lis<string>?), by default None: Either None, or a list of strings with the same
          size as the width of the mutation_matrix representation as a np.ndarray. The strings follow the same
          conventions as the ones in cell_labels.

        Return Values:
        - LabeledMutationMatrix: An object representing mutation_matrix with each row labeled with a label in
          cell_labels and each column labeled with a label in mutation_labels. If either or both the lists of
          labels weren't specified as parameters, the object will be built as if each row/column (or both, 
          depending on the omitted label list) were labeled with a one-based indexing; labels will be strings
          either way.
          The objects fed in as parameters are completely independent from the internal object state, so 
          subsequent modifications to these objects won't alter the content of the LabeledMutationMatrix.
        """

        def is_collection(obj):
            try:
                list(obj)
                return True
            except TypeError:
                return False
        
        # Given how commonly used numpy is, we let its exceptions go through.
        
        mutation_matrix = np.array(mutation_matrix, dtype = 'int')
        if len(mutation_matrix.shape) != 2:
            raise ValueError(
                f"the input matrix should be two-dimensional, but {mutation_matrix} is {len(mutation_matrix.shape)}-dimensional instead."
            )

        if len(mutation_matrix.flat) == 0:
            raise EmptyMatrixError()

        allowed_values = {1, 2, 0}
        for element in mutation_matrix.flat:
            if element not in allowed_values:
                raise ValueError(element)
        height, width = tuple(mutation_matrix.shape)
        
        # Automatically make the labels if none are specified.
        if cell_labels is None:
            cell_labels = [str(i) for i in range(1, height + 1)]
        if mutation_labels is None:
            mutation_labels = [str(i) for i in range(1, width + 1)]

        def _validate_labels(label_list, expected_length):
            if not is_collection(label_list):
                raise TypeError('mutation_labels and cell_labels must be collections, but {label_list} is not')
            if any([isinstance(label_list, str_like_type) for str_like_type in {str, bytes}]):
                raise TypeError(f'strings like {label_list} are not accepted as label lists')
            
            if len(label_list) != expected_length:
                raise MatrixLabelSizeMismatch(f'expected {expected_length} labels, but found {len(label_list)} instead.')
            
            for label in label_list:
                if not isinstance(label, str):
                    raise TypeError(f'each label must be a string, but {label} is not.')
                if re.match('\s', label):
                    raise ValueError(f'labels cannot contain whitespace, newline or other tabulation characters, but {label} does.')

            if len(set(label_list)) != len(label_list):
                raise DuplicateLabelsError()

        _validate_labels(cell_labels, height)
        _validate_labels(mutation_labels, width)
        
        self._data = DataFrame(mutation_matrix, index = cell_labels, columns = mutation_labels, dtype = int, copy = True)

    # The initial choice is to make the matrix immutable and to use external representations that are as general
    # as possible.
        
    def matrix(self):
        """
        Returns a copy of the matrix that was used to initialize the object as a list of lists, where each list is a row 
        in the matrix.

        Parameters: none.

        Returns: 
        - list(list(int)): A list that represents the matrix. The coefficient in position [i][j] will be a 0 if
          the i-th cell doesn't exhibit the j-th mutation, 1 if it does, and 2 if no conclusion could be drawn in the
          analysis. The returned object can be altered without altering the method owner's state.
        """
        return self._data.to_numpy().tolist()
    
    @property
    def cell_labels(self):
        """
        Returns an ordered list of the labels for the cells in the matrix.
        """
        return list(self._data.index)
    
    @property
    def mutation_labels(self):
        """
        Returns an ordered list of the labels for the mutations in the matrix.
        """
        return list(self._data.columns)

    def cells(self):
        """
        Returns a mapping between each cell label and a copy of the corresponding row of the matrix. 

        Returns:
        - dict: a dictionary where each key is a cell label and each value is a copy of the row of 
          the matrix that represents the informations about each mutation for the cell with that label;
          the coefficients of the list are in the same order as the labels in mutation_labels. 
          The list can be intended as a point in the space of the mutations.
        """
        return {lb : list(self._data.loc[lb, :]) for lb in self.cell_labels}

    def mutations(self):
        """
        Returns a mapping between each mutation label and a copy of the corresponding column in the label.

        Parameters: none.

        Returns:
        - dict: a dictionary where each key is a mutation label and each value is a copy of the column of
          the matrix that represents the information for the mutation in each cell of the sample; the coefficients
          of the list are in the same order as the labels in  cell_labels.
        """
        return {lb : list(self._data.loc[:, lb]) for lb in self.mutation_labels}

    # Representing the serializable form with strings lets us use existing formats for the representations.
    def to_serializable_dict(self):
        """
        Dumps the matrix into a dictionary representation which can be easily serialized into a YAML document.
        
        Parameters: none.

        Returns: 
        - dict: a dictionary with the following contents:
          * mutation_matrix: the representation of the matrix as a SASC format string, in which each row is
            separated by a newline and each cell is separated by a space.
          * cell_labels: the representation of the cell labels list. Labels are separated with newlines.
          * mutation_labels: the representation of the mutation labels list. Labels are separated with newlines.
        """
        mutation_matrix_as_string = (
            '\n' +'\n'.join([' '.join([str(el) for el in line]) for line in self.matrix()])
        )
        out = {
            'mutation_matrix' : mutation_matrix_as_string,
            'cell_labels' : '\n' + '\n'.join(self.cell_labels),
            'mutation_labels' : '\n' + '\n'.join(self.mutation_labels)
        }
        return out

    def _from_strings(mutation_matrix, cell_labels = None, mutation_labels = None, matstring_format = "SASC", **kwargs):
        """
        Builds a LabeledMutationMatrix from the string representation of its components, if the representations
        are valid.

        Parameters:
        - mutation_matrix (string): represents a mutation matrix in SASC, SCITE or SPHYR format. 
        - cell_labels (string?), by default None: represents the labels for the cells as a string
          in which the labels are separated by newlines. 
        - mutation_labels (string?): represents the labels for the mutation with the same convention.
        - matstring_format (string in {"SASC", "SPHYR", "SCITE"}): the specification for the format of the
          string.

        Returns:
        - LabeledMutationMatrix: the object representation of the input; refer to the initializer for the
          validation conditions and the behaviour if label files are omitted.
        """
        def _parse_labels(labels_string):
            if labels_string is None:
                return None
            return [lb.strip() for lb in labels_string.splitlines() if lb.strip() != ""]
        
        return LabeledMutationMatrix(
            MatrixParser._from_default_format(matstring_format).parse_matrix(mutation_matrix),
            _parse_labels(cell_labels),
            _parse_labels(mutation_labels)
        )

    def from_serializable_dict(dict_representation):
        """
        Builds a LabeledMutationMatrix from its dict form as it'd be obtained from to_serializable_dict.
        """
        return LabeledMutationMatrix._from_strings(matstring_format = 'SASC', **dict_representation)

    def from_files(matrix_file, matstring_format = 'SASC', cells_file = None, mutations_file = None):
        """
        Reads a matrix file and (facultatively) label files, then uses their content as strings to build
        a LabeledMutationMatrix with the same behaviour as read_from_strings.
        """
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

    def to_files(self, matrix_file, cells_file = None, mutations_file = None):
        """
        Dumps a matrix to a file with a specified format. Also dumps the cell labels and the mutation labels
        if outputs are specified for them as well.
        """

        output_files = [file for file in {matrix_file, cells_file, mutations_file} if file is not None]
        if len(set(output_files)) != len(output_files):
            raise AttributeError('the same file was specified for more than one output.')

        matrix_dict = self.to_serializable_dict()
        with open(matrix_file, 'w+') as f:
            f.write(matrix_dict['mutation_matrix'])
        if cells_file is not None:
            with open(cells_file, 'w+') as f:
                f.write(matrix_dict['cell_labels'])
        if mutations_file is not None:
            with open(mutations_file, 'w+') as f:
                f.write(matrix_dict['mutation_labels'])


class MatrixFileFormat:
    """
    Specifications for the file formats in which matrixes can be stored.

    Attributes:
    - grammar (string): represents a PEG in TatSu format, in which the rule
      labeled with 'matrix' represents a matrix, the rule labeled with 'row' represents
      a matrix row/column (if transpose is true), and the rule 'mcell' represents
      a coefficient in the matrix. The grammar is not validated, so adding new grammars
      must be followed by torough testing to ensure it works as intended.
    - value_map (dict): a map that associates each possible production of the 'cell' rule
      to the value in SASC notation that represents the same information about a mutation; 
      0 means the mutation hasn't occurred for a given cell, 1 means it did and 2 means there's
      not enough evidence to choose.
    - transpose (boolean): if true, then the matrix that is obtained by parsing the file
      will have to be transposed before being fed to the LabeledMutationMatrix initializer.
    """
    def __init__(self, grammar, value_map, transpose = False):
        self.grammar = grammar
        self.value_map = value_map
        self.transpose = transpose

# Should this go into its own file?
_formats = {
    'SASC' : MatrixFileFormat(
        grammar = (
            r"""@@grammar :: SASC
            @@whitespace ::  /(?s)[ \t\r\f\v]+/
            start = matrix $ ;
            mcell = "0" ~ | "1" ~ | "2" ~ ;
            row = { mcell };
            matrix = ("\n").{ row } ;
            """
        ),
        value_map = {'0' : 0, '1': 1, '2' : 2},
        transpose = False
    ),
    
    'SCITE' : MatrixFileFormat(
        grammar = (
            r"""
            @@grammar :: SCITE
            @@whitespace ::  /(?s)[ \t\r\f\v]+/
            start = matrix $ ;
            mcell = "0" ~ | "1" ~ | "2" ~ | "3" ~ ;
            row = { mcell };
            matrix = ("\n").{ row } ;
            """
        ),
        value_map = {'0' : 0, '1' : 1, '2' : 2, '3' : 2},
        transpose = False
    ),
    
    'SPHYR' : MatrixFileFormat(
        grammar = (
            r"""
            @@grammar :: SPHYR
            @@whitespace ::  /(?s)[ \t\r\f\v]+/
            @@eol_comments :: /#([^\n]*?)$/
            start = cell_number_line "\n" snv_number_line "\n" matrix $ ;
            cell_number_line = /\d+/;
            snv_number_line = /\d+/;
            mcell = "0" ~ | "1" ~ | "-1" ~ ;
            row = { mcell };
            matrix = ("\n").{ row };
            """
        ),
        value_map = {'0' : 0, '1': 1, '-1' : 2},
        transpose = True
    )
}


class NotAMatrixError(Exception): pass
class MeaninglessCellError(Exception): pass

class _MatrixSemantics:
    """
    Semantics for a matrix file format. This lets the parsing work as usual
    while also allowing to build a list of each matrix that was found in the string;
    the semantics rely on the names for the rules as specified in the documentation
    for MatrixFileFormat, so new grammars will have to follow them as well.
    """
    def __init__(self, value_map, transpose):
        self.matrix_list = []
        self.current_matrix_rows = []
        self.current_row = []
        self._value_map = value_map.copy()
        self._transpose = transpose
    
    def matrix(self, ast):
        matrix_width = len(self.current_matrix_rows[0])
        if any([len(row) != matrix_width for row in self.current_matrix_rows]):
            raise NotAMatrixError()

        matrix = np.array(self.current_matrix_rows)
        self.current_matrix_rows = []
        if self._transpose: mutation_matrix = np.transpose(matrix)
        self.matrix_list.append(matrix)
        return ast
    
    def row(self, ast):
        if self.current_row != []:
            self.current_matrix_rows.append(self.current_row.copy())
        self.current_row = []
        return ast
    
    def mcell(self, ast):
        try:
            self.current_row.append(self._value_map[ast])
        except KeyError:
            raise MeaninglessCellError(ast)
        return ast

    def default(self, ast):
        return ast

class MatrixParser:
    """
    A parser that reads a matrix string with the specified format and outputs the first
    matrix that is found in the parsed string.
    """
    def __init__(self, file_format):
        self.inner_parser = ts.compile(file_format.grammar)
        self.matrix_builder = _MatrixSemantics(file_format.value_map, file_format.transpose)

    def parse_matrix(self, matrix_string):
        
        self.inner_parser.parse(matrix_string, semantics = self.matrix_builder)

        matrix = self.matrix_builder.matrix_list[0]
        return matrix

    def _from_default_format(matstring_format):
        return MatrixParser(_formats[matstring_format])