import setuptools as st
from Cython.Build import cythonize

sascdir = 'submodules/sasc/'

st.setup(
    name = 'phylopipelinetestbuild',
    ext_modules = cythonize(
        st.Extension("sasc",
                  [sascdir + filepath for filepath in
                      ["bindings/sascpy.pyx",
                       "bindings/sasc-compute.c",
                       "sasc/tree.c",
                       "sasc/mt19937ar.c",
                       "sasc/vector.c",
                       "sasc/utils.c",
                       "sasc/sastep.c"]],
                  define_macros = [],
                  include_dirs = [sascdir + "bindings", sascdir + "sasc"],
                  extra_compile_args = ["-O3", "-fopenmp", "-DNDEBUG"],
                  extra_link_args = ["-O3", "-fopenmp", "-DNDEBUG"]
        ),
        language_level = 3
    ),
    packages = st.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
    ],
    install_requires = [
        'pandas',
        'tatsu'
    ]
)
