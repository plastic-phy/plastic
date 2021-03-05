from setuptools import Extension, setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize(
        Extension("sasc",
                  ["bindings/sascpy.pyx",
                   "bindings/sasc-compute.c",
                   "sasc/tree.c",
                   "sasc/mt19937ar.c",
                   "sasc/vector.c",
                   "sasc/utils.c",
                   "sasc/sastep.c"],
                  define_macros = [("DNDEBUG",)],
                  include_dirs = ["bindings", "sasc"],
                  extra_compile_args = ["-O3", "-fopenmp"],
                  extra_link_args = ["-fopenmp"]
        ),
        language_level = 3
    )
)
