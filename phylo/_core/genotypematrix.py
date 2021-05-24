from pandas import DataFrame
import numpy as np
import tatsu as ts
import re
from copy import deepcopy


class DuplicateLabelsError(Exception): pass


class EmptyMatrixError(Exception): pass


class MatrixLabelSizeMismatch(Exception): pass


class GenotypeMatrix:

    def __init__(self, genotype_matrix, cell_labels=None, mutation_labels=None):
        """
        Validates a matrix that describes the results of an analysis of the mutations that have
        occurred in a sample of cells.

        Parameters:
            genotype_matrix:
                Any twice-subscriptable collection of values that can be converted to a
                bidimensional, nonempty np.ndarray of integer in which each value is either 0, 1 or 2.
                Each field in the matrix represents information on whether if a mutation has or hasn't occurred
                in a certain cell. 0 = no mutation, 1 = mutation, 2 = unknown.
            cell_labels (list<string>?), by default None:
                Either None, or a list of strings with the same size as the height of the genotype matrix.
                If the list exists, each element must be a non-empty string without spaces, newlines
                or other tabulation characters.
                If a list is provided, then the cells will be labeled in that order. Otherwise, an incremental
                labeling scheme will be used.
            mutation_labels (list<string>?), by default None:
                Either None, or a list of strings with the same size as the width of the genotype matrix.
                If the list exists, each element must be a non-empty string without spaces, newlines
                or other tabulation characters.
                If a list is provided, then the mutations will be labeled in that order. Otherwise, an incremental
                labeling scheme will be used.

        Returns:
            GenotypeMatrix:
                An object representing the matrix, labeled as specified above. The object's internal state
                is completely independent from the objects used for its initialization.
        """

        def is_collection(obj):
            try:
                list(obj)
                return True
            except TypeError:
                return False

        # Given how commonly used numpy is, we let its exceptions go through.
        genotype_matrix = np.array(genotype_matrix, dtype='int')
        if len(genotype_matrix.shape) != 2:
            raise ValueError(
                f"the input matrix should be two-dimensional, but {genotype_matrix} is {len(genotype_matrix.shape)}-dimensional instead."
            )

        if len(genotype_matrix.flat) == 0:
            raise EmptyMatrixError()

        allowed_values = {1, 2, 0}
        for element in genotype_matrix.flat:
            if element not in allowed_values:
                raise ValueError(element)
        height, width = tuple(genotype_matrix.shape)

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
                raise MatrixLabelSizeMismatch(
                    f'expected {expected_length} labels, but found {len(label_list)} instead.')

            for label in label_list:
                if not isinstance(label, str):
                    raise TypeError(f'each label must be a string, but {label} is not.')
                if re.search(r'\s', label):
                    raise ValueError(
                        f'labels cannot contain whitespace, newline or other tabulation characters, but {label} does.')

            if len(set(label_list)) != len(label_list):
                raise DuplicateLabelsError()

        _validate_labels(cell_labels, height)
        _validate_labels(mutation_labels, width)

        self._data = DataFrame(genotype_matrix, index=cell_labels, columns=mutation_labels, dtype=int, copy=True)

    def matrix(self):
        """
        Returns a copy of the matrix that was used to initialize the object as a list of lists, where each list is a row 
        in the matrix.

        Parameters:
            none.

        Returns:
            np.ndarray:
                A numpy array that represents the matrix. The coefficient in position [i][j] will be a 0 if
                the i-th cell doesn't exhibit the j-th mutation, 1 if it does, and 2 if no conclusion could be drawn
                in the analysis. The returned object is independent from the internal representation.
        """
        return deepcopy(self._data.to_numpy())

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

        Parameters:
            none.

        Returns:
            dict:
                A dictionary where each key is a cell label and each value is a copy of the row of
                the matrix that represents the information about each mutation for the cell with that label;
                the coefficients of the list are in the same order as the labels in mutation_labels.
        """
        return {lb: list(self._data.loc[lb, :]) for lb in self.cell_labels}

    def mutations(self):
        """
        Returns a mapping between each mutation label and a copy of the corresponding column in the label.

        Parameters:
            none.

        Returns:
            dict:
                A dictionary where each key is a mutation label and each value is a copy of the column of
                the matrix that represents the information for the mutation in each cell of the sample;
                the coefficients of the list are in the same order as the labels in  cell_labels.
        """
        return {lb: list(self._data.loc[:, lb]) for lb in self.mutation_labels}

    def with_automatic_mutation_labels(self):
        """
        Returns a copy of the matrix where the mutation labels are replaced
        by a range from 1 to the number of mutations. Useful after mutation
        clusterization if there is a need to shorten the labels and/or to
        ignore their cluster structure.
        """
        return self.with_mutation_labels(None)

    def with_mutation_labels(self, new_mutation_labels):
        """
        Returns a copy of the matrix with a new set of mutation labels.
        """
        return GenotypeMatrix(self.matrix(), cell_labels=self.cell_labels, mutation_labels=new_mutation_labels)

    def with_cell_labels(self, new_cell_labels):
        """
        Returns the same matrix with a new set of cell labels
        """
        return GenotypeMatrix(self.matrix(), cell_labels=new_cell_labels, mutation_labels=self.mutation_labels)

    def to_serializable_dict(self):
        """
        Dumps the matrix into a dictionary-of-strings representation. The matrix is represented in SASC format,
        while the label lists are represented as newline-separated string.
        """
        genotype_matrix_as_string = (
            '\n'.join([' '.join([str(el) for el in line]) for line in self.matrix()])
        )
        out = {
            'genotype_matrix': genotype_matrix_as_string,
            'cell_labels': '\n'.join(self.cell_labels),
            'mutation_labels': '\n'.join(self.mutation_labels)
        }
        return out

    @classmethod
    def _from_strings(cls, genotype_matrix, cell_labels=None, mutation_labels=None, matstring_format="SASC"):
        """
        Builds a GenotypeMatrix from the string representation of its components, if the representations
        are valid.

        Parameters:
            genotype_matrix(string):
                Represents a genotype matrix in SASC, SCITE or SPHYR format.
            cell_labels(string?), by default None:
                Represents the labels for the cells as a string in which the labels are separated by newlines.
                If the string isn't specified the labels will be generated as described in the __init__.
            mutation_labels(string?), by default None:
                Represents the labels for the mutations as a string in which the labels are separated by newlines.
                If the string isn't specified the labels will be generated as described in the __init__.
            matstring_format(string in {"SASC", "SPHYR", "SCITE"}):
                The specification for the format of the matrix string.

        Returns:
            GenotypeMatrix:
                The object representation of the input; refer to the initializer for the
                validation conditions and the behaviour if label files are omitted.
        """

        def _parse_labels(labels_string):
            if labels_string is None:
                return None
            return [lb.strip() for lb in labels_string.splitlines() if lb.strip() != ""]

        return cls(
            _matrix_parser_from_default_format(matstring_format).parse_matrix(genotype_matrix),
            _parse_labels(cell_labels),
            _parse_labels(mutation_labels)
        )

    @classmethod
    def from_serializable_dict(cls, dict_representation):
        """
        Builds a GenotypeMatrix from its dict form as it'd be obtained from to_serializable_dict.
        """
        return cls._from_strings(matstring_format='SASC', **dict_representation)

    @classmethod
    def from_files(cls, matrix_file, matstring_format='SASC', cells_file=None, mutations_file=None):
        """
        Reads a matrix file and (facultatively) label files, then uses their content as strings to build
        a GenotypeMatrix with the same behaviour as read_from_strings.
        """

        def _read_nullable(file_name, default_output=None):
            if file_name is None:
                return default_output
            else:
                with open(file_name, 'r') as f:
                    return f.read()

        with open(matrix_file, 'r') as f:
            genotype_matrix = f.read()

        return cls._from_strings(
            genotype_matrix=genotype_matrix,
            cell_labels=_read_nullable(cells_file),
            mutation_labels=_read_nullable(mutations_file),
            matstring_format=matstring_format
        )

    def to_files(self, matrix_file, cells_file=None, mutations_file=None):
        """
        Dumps a matrix to a file with a specified format. Also dumps the cell labels and the mutation labels
        if outputs are specified for them as well. The files will be created if they don't exist,
        but only if the target directory already exists. If the files already exist, **they will be overwritten**.
        """

        output_files = [file for file in {matrix_file, cells_file, mutations_file} if file is not None]
        if len(set(output_files)) != len(output_files):
            raise AttributeError('the same file was specified for more than one output.')

        matrix_dict = self.to_serializable_dict()
        with open(matrix_file, 'w') as f:
            f.write(matrix_dict['genotype_matrix'])
        if cells_file is not None:
            with open(cells_file, 'w') as f:
                f.write(matrix_dict['cell_labels'])
        if mutations_file is not None:
            with open(mutations_file, 'w') as f:
                f.write(matrix_dict['mutation_labels'])


class MatrixFileFormat:
    """
    Specifications for the file formats in which matrixes can be stored.
    """

    def __init__(self, grammar, value_map, transpose=False):
        """
        Parameters:
            grammar(string):
                Represents a PEG in TatSu format. Some rules have semantics associated to them: the rule
                labeled with 'matrix' represents a matrix, the rule labeled with 'row' represents
                a matrix row/column (depending on the value for transpose), and the rule 'mcell' represents
                a coefficient in the matrix. The grammar is not validated, so adding new grammars
                must be followed by thorough testing to ensure it works as intended.
            value_map(dict):
                A map that associates each possible production of the 'mcell' rule
                to the value in SASC notation that represents the same information about a mutation;
                The usual convention is used: 0 = no mutation, 1 = mutation, 2 = unknown.
            transpose(boolean):
                If true, then the matrix that is obtained by parsing the file
                will have to be transposed before being fed to the GenotypeMatrix initializer.
        """
        self.grammar = grammar
        self.value_map = value_map
        self.transpose = transpose


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
        if self._transpose: matrix = np.transpose(matrix)
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
        self.inner_parser.parse(matrix_string, semantics=self.matrix_builder)

        matrix = self.matrix_builder.matrix_list[0]
        return matrix


_formats = {
    'SASC': MatrixFileFormat(
        grammar=(
            r"""@@grammar :: SASC
            @@whitespace ::  /(?s)[ \t\r\f\v]+/
            start = matrix $ ;
            mcell = "0" ~ | "1" ~ | "2" ~ ;
            row = { mcell };
            matrix = ("\n").{ row } ;
            """
        ),
        value_map={'0': 0, '1': 1, '2': 2},
        transpose=False
    ),

    'SCITE': MatrixFileFormat(
        grammar=(
            r"""
            @@grammar :: SCITE
            @@whitespace ::  /(?s)[ \t\r\f\v]+/
            start = matrix $ ;
            mcell = "0" ~ | "1" ~ | "2" ~ | "3" ~ ;
            row = { mcell };
            matrix = ("\n").{ row } ;
            """
        ),
        value_map={'0': 0, '1': 1, '2': 2, '3': 2},
        transpose=True
    ),

    'SPHYR': MatrixFileFormat(
        grammar=(
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
        value_map={'0': 0, '1': 1, '-1': 2},
        transpose=False
    )
}


def _matrix_parser_from_default_format(matstring_format):
    return MatrixParser(_formats[matstring_format])
