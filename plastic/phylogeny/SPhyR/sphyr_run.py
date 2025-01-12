from plastic import GenotypeMatrix, PhylogenyTree
import subprocess
import os
import tempfile 
from typing import Union

def M_T_check(par):
    if isinstance(par, int) and par == -1 or par > 0:
        return True
    else: 
        return False

def N_check(par):
    if isinstance(par, int) and par > 0:
        return True
    else:
        return False

def k_s_check(par):
    if isinstance(par, int) and par >= 0:
        return True
    else: 
        return False

def a_b_check(par):
    if isinstance(par, float) and par > 0 and par < 1:
        return True
    else: 
        return False

def lC_lT_t_check(par):
    if isinstance(par, int) and par > 0:
        return True
    else: 
        return False

def output_check(par):
    if isinstance(par, str):
        return True
    else:
        return False
    
def path_check(path):
    if isinstance(path, str) and len(path) > 0:
        return True
    else:
        return False

def create_kDPFC_command(kDPFC_path, visualize_path, input_file, cells_file, mutations_file, M, N, T, a, b, k, lC, lT, s, t, output):
    if not path_check(kDPFC_path):
        os.remove(input_file)
        os.remove(cells_file)
        os.remove(mutations_file)
        raise Exception('wrong value for the path to kDPFC:' + str(kDPFC_path))
    if not path_check(visualize_path):
        os.remove(input_file)
        os.remove(cells_file)
        os.remove(mutations_file)
        raise Exception('wrong value for the path to visualize:' + str(visualize_path))
    if not M_T_check(M):
        os.remove(input_file)
        os.remove(cells_file)
        os.remove(mutations_file)
        raise Exception('wrong value for the parameter M:' + str(M))
    if not M_T_check(T):
        os.remove(input_file)
        os.remove(cells_file)
        os.remove(mutations_file)
        raise Exception('wrong value for the parameter T:' + str(T))
    if not N_check(N):
        os.remove(input_file)
        os.remove(cells_file)
        os.remove(mutations_file)
        raise Exception('wrong value for the parameter N:' + str(N))
    if not k_s_check(k):
        os.remove(input_file)
        os.remove(cells_file)
        os.remove(mutations_file)
        raise Exception('wrong value for the parameter k:' + str(k))
    if not k_s_check(s):
        os.remove(input_file)
        os.remove(cells_file)
        os.remove(mutations_file)
        raise Exception('wrong value for the parameter s:' + str(s))
    if not a_b_check(a):
        os.remove(input_file)
        os.remove(cells_file)
        os.remove(mutations_file)
        raise Exception('wrong value for the parameter a:' + str(a))
    if not a_b_check(b):
        os.remove(input_file)
        os.remove(cells_file)
        os.remove(mutations_file)
        raise Exception('wrong value for the parameter b:' + str(b))
    if not lC_lT_t_check(lC):
        os.remove(input_file)
        os.remove(cells_file)
        os.remove(mutations_file)
        raise Exception('wrong value for the parameter lC:' + str(lC))
    if not lC_lT_t_check(lT):
        os.remove(input_file)
        os.remove(cells_file)
        os.remove(mutations_file)
        raise Exception('wrong value for the parameter lT:' + str(lT))
    if not lC_lT_t_check(t):
        os.remove(input_file)
        os.remove(cells_file)
        os.remove(mutations_file)
        raise Exception('wrong value for the parameter t:' + str(t))
    if not output_check(output):
        os.remove(input_file)
        os.remove(cells_file)
        os.remove(mutations_file)
        raise Exception('the name of the output file cannot be empty:' + str(output))
    
    command = [kDPFC_path, '-M', f'{M}', '-N', f'{N}', '-T', f'{T}', '-a', f'{a}', '-b', f'{b}',
            '-k', f'{k}', '-lC', f'{lC}', '-lT', f'{lT}', '-s', f'{s}', '-t', f'{t}', input_file]
    
    if len(output) == 0:
        command.append('output')
    else:
        command.append(output)

    return command

def run_sphyr(matrix : GenotypeMatrix, kDPFC_path : str, visualize_path : str, wd : str = os.getcwd(), M : int = -1, N : int = 10,
              T : int = -1, a : float = 1e-3, b : float = 0.3, k : int = 1, lC : int = 15, lT : int = 10, s : int = 0,
              t : int = 1, name_kDPFC : str = 'output', matrix_out : bool = False) -> PhylogenyTree | list[PhylogenyTree, GenotypeMatrix]:
    """
    Returns a PhylogenyTree object or a list containing the PhylogenyTree object in position 0 and a GenotypeMatrix in position 1:
    The PhylogenyTree object is the result of the kDPFC and visualize programs
    The GenotypeMatrix object is the translation of the output of kDPFC program
    NOTE that SPhyR must be downloaded and compiled by the user.
    For better informations read the README file for SPhyR.

    input:
        matrix: GenotypeMatrix object used to run kDPFC
        kDPFC_path: path to the kDPFC executable file included
        visualize_path: path to the visualize executable file included
        wd: working directory set by default to os.getcwd()
        M: int representing the memory limit in MB, by default set to -1
        N: int representing the number of restarts, by default set to 10
        T: int representing the time limit in seconds, by default set to -1 
        a: float representing the false positive rate, by default set to 1e-3
        b: float representing false negative rate, by default set to 0.3
        k: int representing the maximum number of losses per SNV, by default set to 1
        lC: int representing the number of character clusters, by default set to 15
        lT: int representing the number of taxon clusters, by default set to 10
        s: int representing the random number generator seed, by default set to 0
        t: int representing the number of threads, by default set to 1
        name_kDPFC: string representing the name of the output of the kDPFC program, by default set to 'output'
        matrix_out: boolean representing if the matrix has to be retuned or not, if yes the list described above is returned, by default set to False
    """

    if len(name_kDPFC) == 0:
        name_kDPFC = 'output'

    with tempfile.NamedTemporaryFile(mode='w+', dir=f'{wd}', delete=False) as matrix_file, \
         tempfile.NamedTemporaryFile(mode='w+', dir=f'{wd}', delete=False) as cells_file, \
         tempfile.NamedTemporaryFile(mode='w+', dir=f'{wd}', delete=False) as mutations_file:
         
        matrix_file_path = matrix_file.name
        cells_file_path = cells_file.name
        mutations_file_path = mutations_file.name

        matrix.to_files(matrix_file = matrix_file_path, cells_file = cells_file_path, mutations_file = mutations_file_path,
                        output_type = 'SPhyR')

        kDPFC_command = create_kDPFC_command(kDPFC_path=kDPFC_path, visualize_path=visualize_path, input_file=matrix_file_path, cells_file=cells_file_path, 
                                             mutations_file=mutations_file_path, M=M, N=N, T=T, a=a, b=b, k=k, lC=lC, lT=lT,
                                             s=s, t=t, output=name_kDPFC)
        
        print('starting kDPFC process...')
        
        process = subprocess.Popen(kDPFC_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = f'{wd}')
        stdout, stderr = process.communicate()

        if process.returncode == 0:
            output_list = [None, None]
            if matrix_out:
                output_list[1] = GenotypeMatrix.from_files(matrix_file = f'{wd}/output', cells_file = cells_file_path, mutations_file=mutations_file_path, str_parser = 'KDPFCParser')

            os.remove(matrix_file_path)

            visualize_command = [visualize_path, '-c', mutations_file_path,'-t', cells_file_path, name_kDPFC]
            print('starting visualize process...')

            process = subprocess.Popen(visualize_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = f'{wd}')
            stdout, stderr = process.communicate()

            os.remove(cells_file_path)
            os.remove(mutations_file_path)

            if process.returncode == 0:
                print('creating the tree...')
                tree = PhylogenyTree.from_dotstring(stdout.decode())
                
                if matrix_out:
                    output_list[0] = tree
                    return output_list
                else:
                    return tree
            
            else:
                raise Exception("Error during the SPhyR visualize execution, error number: " + str(process.returncode) + 
                                ", error message: " + stderr.decode())

        else:
            os.remove(matrix_file_path)
            os.remove(cells_file_path)
            os.remove(mutations_file_path)
            raise Exception("Error during the SPhyR kDPFC execution, error number: " + str(process.returncode) + 
                            ", error message: " + stderr.decode())
