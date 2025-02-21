from plastic import GenotypeMatrix, PhylogenyTree
import os
import tempfile
import subprocess

def search_matrix_placeholder(matrix_placeholder, parameters):
    if isinstance(matrix_placeholder, str) and len(matrix_placeholder) > 0:
        presence = False
        for command in parameters:
            if matrix_placeholder in command:
                presence = True
                break
        if not presence:
            raise Exception('there is not the specified matrix place holder in the list of parameters')
    else:
        raise Exception('the specified matrix place holder is not accepted')

def cells_mutations_placeholders_check(matrix_placeholder, cells_placeholder, mutations_placeholder, parameters):
    search_cells = False
    search_mutations = False
    if cells_placeholder is not None:
        if isinstance(cells_placeholder, str) and len(cells_placeholder) > 0:
            if cells_placeholder != matrix_placeholder:
                search_cells = True
            else:
                raise Exception('cells place holder has the same name of matrix place holder')
        else:
            raise Exception('cells place holder is not accepted')

    if mutations_placeholder is not None:
        if isinstance(mutations_placeholder, str) and len(mutations_placeholder) > 0:
            if mutations_placeholder != matrix_placeholder:
                if search_cells:
                    if mutations_placeholder != cells_placeholder:
                        search_mutations = True
                    else:
                        raise Exception('mutations place holder has the same name of cells place holder')
            else:
                raise Exception('mutations place holder has the same name of matrix place holder')
        else:
            raise Exception('mutations place holder is not accepted')
    
    for command in parameters:
        if not search_cells and not search_mutations:
            break
        if search_cells and cells_placeholder in command:
            search_cells = False
        if search_mutations and mutations_placeholder in command:
            search_mutations = False

    if search_cells or search_mutations:
        raise Exception('error: at least one between the cells place holder or the mutations place holder cannot be found in the list of command')

def create_command(parameters, matrix, cells, mutations, matrix_file_path, cells_file_path,
                   mutations_file_path) -> list[list[str]]:
    final_commands = []
    for command in parameters:
        internal_command = []
        for e in command:
            if e == matrix:
                internal_command.append(matrix_file_path)
            elif cells is not None and e == cells:
                internal_command.append(cells_file_path)
            elif mutations is not None and e == mutations:
                internal_command.append(mutations_file_path)
            else:
                internal_command.append(e)
        final_commands.append(internal_command)
    return final_commands

def inference(matrix : GenotypeMatrix, parameters : list[list[str]]
                      , matrix_placeholder : str,
                      output : str, wd : str = os.getcwd(), cells_placeholder : str = None,
                      mutations_placeholder : str = None, type_tree = PhylogenyTree) -> PhylogenyTree:
    
    '''
    Returns a PhylogenyTree object as the result of executing a sequence of commands specified in parameters, all the executables must be 
    downloaded and compiled by the user before the use.

    input:
        matrix: GenotypeMatrix or a subclass object that implements a correct version of the to_files method
        parameters: list of lists of strings (list[list[str]]) where the strings are all the parts of the commands and the internal lists are 
            the total commands, in at least one of the internal lists must be the value of matrix_placeholder
            example:
                if the command to run the executable is:
                path/to/the/executable/included path/to/the/matrix/file -option1 a -option2 b -option3 c
                than the parameters variable must be:
                parameters = [['path/to/the/executable/included', 'matrix_placeholder', '-option1', 'a', '-option2', 'b', '-option3', 'c']]
                if two or more different commands are needed than the parameters variable must be:
                parameters = [['path/to/the/first/executable/included', 'matrix_placeholder_if_needed_in_this_command', '-option1', 'a', '-option2', 'b', '-option3', 'c'],
                                'path/to/the/second/executable/included', 'matrix_placeholder_if_needed_in_this_command', '-option1', 'd', '-option2', 'e', '-option3', 'f']
        matrix_placeholder: str representing the placeholder for the matrix in parameters, this value will be replaced before the use
            with the path to the temporary file used to store the matrix data, must be unique within parameters to prevent errors
        output: str representing the name of the output obtained as result of the execution of the last command to be 
            used to construct the PhylogenyTree. if the output is printed to the terminal use 'stdout'
        wd: str representing the path to the working directory
        cells_placeholder: str representing the placeholder for the cells file in parameters, this value will be replaced before the use
            with the path to the temporary file used to store the cells names, must be unique and different from matrix_placeholder
        mutations_placeholder: str representing the placeholder for the mutations file in parameters, this value will be replaced before the use
            with the path to the temporary file used to store the mutations' names, must be unique and different from both matrix_placeholder and cells_placeholder
        type_tree: PhylogenyTree or a subclass implementing the from_file or from_dotstring method used to instantiate the tree that will
            be returned, type_tree must not be an object of a class, it must be the class
    '''
    
    if not isinstance(output, str) or len(output) < 1:
        raise Exception('the output string is incorrect')
    
    if type_tree is None or not issubclass(type_tree, PhylogenyTree):
        raise Exception('type tree is not valid')
    
    search_matrix_placeholder(matrix_placeholder=matrix_placeholder, parameters=parameters)

    cells_mutations_placeholders_check(matrix_placeholder=matrix_placeholder, cells_placeholder=cells_placeholder,
                                       mutations_placeholder=mutations_placeholder, parameters=parameters)

    with tempfile.NamedTemporaryFile(mode='w+', dir=f'{wd}', delete=False) as matrix_file, \
         tempfile.NamedTemporaryFile(mode='w+', dir=f'{wd}', delete=False) as cells_file, \
         tempfile.NamedTemporaryFile(mode='w+', dir=f'{wd}', delete=False) as mutations_file:

        matrix_file_path = matrix_file.name
        cells_file_path = cells_file.name
        mutations_file_path = mutations_file.name

        matrix.to_files(matrix_file_path, cells_file_path, mutations_file_path)

        final_commands = create_command(parameters=parameters, matrix = matrix_placeholder, cells = cells_placeholder,
                                        mutations = mutations_placeholder, matrix_file_path= matrix_file_path,
                                        cells_file_path= cells_file_path, mutations_file_path= mutations_file_path)

        for command in range(len(final_commands)):
            process = subprocess.Popen(final_commands[command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd = wd)
            stdout, stderr = process.communicate()
            if process.returncode == 0:
                if command == len(final_commands) - 1:
                    if output != 'stdout':
                        tree = type_tree.from_file(f'{wd}/{output}')
                    else:
                        tree = type_tree.from_dotstring(stdout.decode())
                    os.remove(matrix_file_path)
                    os.remove(cells_file_path)
                    os.remove(mutations_file_path)
                    return tree
            else:
                os.remove(matrix_file_path)
                os.remove(cells_file_path)
                os.remove(mutations_file_path)
                print(stderr.decode())
                raise Exception('error during the execution of process number ' + str(command) + ': ' + str(final_commands[command]) +
                                ', number error '+ str(process.returncode) + ' has been reported: ' + stderr.decode())