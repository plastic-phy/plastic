import subprocess
import tempfile
import os
from plastic import GenotypeMatrix, PhylogenyTree

def r_l_check(par):
    """
    checks that par is int and greater than 0
    """
    if isinstance(par, int) and par > 0:
        return True
    else:
        return False

def seed_check(par):
    """
    checks that par is int greater or equal to 0
    """
    if isinstance(par, int) and par >= 0:
        return True
    else:
        return False

def ad_fd_check(par):
    """
    checks that par is float and greater than 0
    """
    if isinstance(par, float) and par > 0:
        return True
    else:
        return False

def s_transpose_a_check(par):
    """
    checks that par is bool
    """
    if isinstance(par, bool):
        return True
    else:
        return False

def e_check(par):
    """
    checks that par is float between 0 and 1 , -1 is default
    """
    if isinstance(par, float) and par >= 0 and par <= 1 or par == -1.0:
        return True
    else:
        return False

def x_sd_check(par):
    """
    checks that par is float greater than 0
    """
    if isinstance(par, float) and par > 0:
        return True
    else: 
        return False

def o_t_check(par):
    """
    checks that par is a not empty string
    """
    if isinstance(par, str):
        return True
    else:
        return False

def g_check(par):
    """
    checks that par is float
    """
    if isinstance(par, float):
        return True
    else:
        return False

def move_probs_check(par, val):
    """
    checks that par is a tuple of 3 floats if val is false, 2 otherwise, all the elements in par must be values between 0 and 1
    """
    if isinstance(par, tuple):
        if val and len(par) == 2:
            for e in par:
                if not isinstance(e, float) or e <= 0 or e >= 1:
                    return False
            return True
        elif not val and len(par) == 3:
            for e in par:
                if not isinstance(e, float) or e <= 0 or e >= 1:
                    return False
            return True
        else: 
            return False
    else:
        return False

def path_check(path):
    if isinstance(path, str) and len(path) > 0:
        return True
    else:
        return False

def insert_parameters(i, scite_path, names, n, m, ad, fd, r, l, s,
                   transpose, e, x, sd, o, a, g, seed, t, move_probs):
    command = [scite_path, '-i', f'{i}', '-n', f'{n}', '-m', f'{m}', '-r', f'{r}',
               '-l', f'{l}', '-fd', f'{fd}', '-ad', f'{ad}']

    if s: 
        command.append('-s')
    if transpose:
        command.append('-transpose')
    if e != -1:
        command.extend(('-e', f'{e}'))
    command.extend(('-x', f'{x}', '-sd', f'{sd}'))
    if len(o) == 0:
        command.extend(('-o', 'output'))
    elif len(o) > 0:
        command.extend(('-o', f'{o}'))
    command.extend(('-names', f'{names}'))
    if a:
        command.append('-a')
    command.extend(('-max_treelist_size', '1', '-g', f'{g}'))
    if seed != 0:
        command.extend(('-seed', f'{seed}'))
    if len(t) > 0:
        command.extend(('-t', f'{t}'))
    command.append('-move_probs')
    for i in move_probs:
        command.insert(len(command), str(i))

    return command
    
def create_command(i, scite_path, names, n, m, fd, ad, r, l, s,
                   transpose, e, x, sd, o, a, g, seed, t, move_probs)-> list:
    if not path_check(scite_path):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the scite_path:' + str(scite_path))
    if not r_l_check(r):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter r:' + str(r))
    if not r_l_check(l):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter l:' + str(l))
    if not seed_check(seed):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter seed:' + str(seed))
    if not ad_fd_check(fd):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter fd:' + str(fd))
    if not ad_fd_check(ad):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter ad:' + str(ad))
    if not s_transpose_a_check(s):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter s:' + str(s))
    if not s_transpose_a_check(transpose):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter transpose:' + str(transpose))
    if not s_transpose_a_check(a):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter a:' + str(a))
    if not e_check(e):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter e:' + str(e))
    if not x_sd_check(x):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter x:' + str(x))
    if not x_sd_check(sd):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter sd:' + str(sd))
    if not o_t_check(o):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter o:' + str(o))
    if not o_t_check(t):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter t:' + str(t))
    if not g_check(g):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter g:' + str(g))
    if not move_probs_check(move_probs, transpose):
        os.remove(i)
        os.remove(names)
        raise Exception('wrong value for the parameter move_probs:' + str(move_probs))
    
    return insert_parameters(i = i, scite_path = scite_path, names = names, n = n, m = m, fd = fd,
                             ad = ad, r = r, l = l, s = s, transpose = transpose, e = e, x = x,
                             sd = sd, o = o, a = a, g = g, seed = seed, t = t, move_probs = move_probs)

def run_scite(matrix : GenotypeMatrix, scite_path : str, fd : float, ad : float, wd : str = os.getcwd(), r : int = 1,
              l : int = 900000, s : bool = False, transpose : bool = False, e : float = -1.0, x : float = 10.0, sd : float = 0.1,
              o : str = 'output', a : bool = False, g : float = 1.0, seed : int = 0, t : str = '',
              move_probs : tuple = (0.55, 0.4, 0.05))-> PhylogenyTree: 
    """
    Returns a PhylogenyTree object as the best result of the downloaded and compiled SCITE program by the user.
    For better informations read the README file for SCITE

    input:
        matrix: GenotypeMatrix object used to run SCITE
        scite_path: path to the SCITE executable file included
        fd (false discoveries): float value representing the estimates false positive rate of the sequencing experiment
        ad (allelic dropout): float value representing the estimated false negative rate of the sequencing experiment
        wd: working directory set by default to os.getcwd()
        r: int representing the desired number of repetition of the MCMC, by default set to 1
        l: int representing the desired chain length of each MCMC repetition, by default set to 900000
        s: bool representing if the option '-s' is selected for the program, this option causes the sample attachment
           points to be marginalized, by default set to False
        transpose: bool representing if the option '-transpose' is selected, this option changes the tree representation
           from mutation tree to rooted binary leaf-labelled tree, by default set to False
        e: float between 0 and 1 representing the learning of error beta, specify the probability to chose the move for
           changing the error rate in the MCMC, by default set to -1.0 to be ignored by the program
        x: float representing the scaling of the known error rate for the MH jump, by default set to 10
        sd: float representing the prior standard deviation for AD erro rate, by default set to 0.1
        o: string representing the name used by SCITE as base name for the output files
        a: bool representing if the option '-a' is selected for the program, with this option SCITE adds the individual 
           cells as additional nodes to the reported trees, read the SCITE README file for more informations
        g: float representing the desired value of gamma, read the SCITE README file for more informations
        seed: int used as a fixed seed for the random number generator
        t: string representing the path to a file with the true tree in GraphViz format
        move_probs: tuple of three floats which changes the default probabilities for the three MCMC moves to the specified 
                    values. If -transpose is selected and move_probs is provided then move_probs must be a tuple of two floats,
                    read the SCITE README file for more informations

        NOTE: the '-p' and '-no_tree_list' options are not available
    """
    if transpose:
        if move_probs == (0.55, 0.4, 0.05):
            move_probs = (0.4, 0.6)

    with tempfile.NamedTemporaryFile(mode='w+', dir=f'{wd}', delete=False) as matrix_file, \
         tempfile.NamedTemporaryFile(mode='w+', dir=f'{wd}',delete=False) as mutations_file:
         
        matrix_file_path = matrix_file.name
        mutations_file_path = mutations_file.name

        matrix.to_files(matrix_file_path, mutations_file = mutations_file_path, output_type = 'scite')

        n_m = (matrix._data.shape[1], matrix._data.shape[0])
        command = create_command(i=matrix_file_path, scite_path=scite_path, names=mutations_file_path, n = n_m[0], 
                                m = n_m[1], fd = fd, ad = ad, r = r, l = l, s = s, transpose = transpose, e = e, 
                                x = x, sd = sd, o = o, a = a, g = g, seed = seed, t = t, move_probs = move_probs)
        print('starting subprocess...')
        process = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = f'{wd}')
        stdout, stderr = process.communicate()

        os.remove(matrix_file_path)
        os.remove(mutations_file_path)
        
        if process.returncode == 0:
            print('creating the tree...')
            tree = PhylogenyTree.from_file(f'{wd}/{o}_ml0.gv')
            return tree
        else:
            raise Exception("Error during the scite execution, error number: " + str(process.returncode) + ", error message: " +
                            stderr.decode())