import setuptools as st
from Cython.Build import cythonize
import numpy

sascdir = 'plastic/phylogeny/sasc/'
with open('README.md', 'r') as f:
    long_description = f.read()

extensions = [
    st.Extension(
        "plastic.phylogeny.sasc._sasc",
        [sascdir + "bindings/_sasc.pyx"] +  # Cython source file
        [sascdir + filepath for filepath in [
            "bindings/sasc-compute.c",
            "sasc/tree.c",
            "sasc/mt19937ar.c",
            "sasc/vector.c",
            "sasc/utils.c",
            "sasc/sastep.c"
        ]],
        include_dirs=[sascdir + "bindings", sascdir + "sasc", numpy.get_include()],
        extra_compile_args=["-O3", "-fopenmp", "-DNDEBUG"],
        extra_link_args=["-O3", "-fopenmp", "-DNDEBUG"]
    )
]

st.setup(
    name="plastic",
    version="0.0.2",
    description="Simpler and Faster Development of Tumor Phylogeny Pipelines",
    author="Lorenzo Lucarella, Simone Ciccolella, Andy Giacon",
    author_email="l.lucarella@campus.unimib.it, simone.ciccolella@unimib.it, a.giacon@campus.unimib.it",
    long_description_content_type="text/markdown",
    long_description=long_description,
    url="https://github.com/plastic-phy/plastic",
    packages=st.find_packages(),
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    install_requires=[
        "pandas",
        "tatsu",
        "mp3treesim",
        "pygraphviz",
        "kmodes",
        "networkx",
        "colour",
        "numpy",
        "matplotlib"
    ],
    python_requires=">=3.10",
    ext_modules=cythonize(extensions, language_level="3"),
)

