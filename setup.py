import setuptools as st
from Cython.Build import cythonize

sascdir = 'phylo/sasc/'
with open('README.md', 'r') as f:
    long_description = f.read()

st.setup(
    name="phylopipeline",
    version="0.0.1",
    description="Pipeline for the AlgoLab suite of cancer phylogeny inference tools",
    author="Lorenzo Lucarella",
    author_email="l.lucarella@campus.unimib.it",
    long_description_content_type="text/markdown",
    long_description=long_description,
    url="https://github.com/AlgoLab/phylo-pipeline/",
    packages=st.find_packages(),
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    install_requires=[
        'pandas',
        'tatsu',
        'mp3treesim',
        'pygraphviz',
        'kmodes',
        'networkx',
        'colour',
        'numpy'
    ],
    python_requires=">=3.8",
    ext_modules=[
        st.Extension(
            "phylo.sasc",
            [sascdir + filepath for filepath in
             ["bindings/sasc.c",
              "bindings/sasc-compute.c",
              "sasc/tree.c",
              "sasc/mt19937ar.c",
              "sasc/vector.c",
              "sasc/utils.c",
              "sasc/sastep.c"]],
            include_dirs=[sascdir + "bindings", sascdir + "sasc"],
            extra_compile_args=["-O3", "-fopenmp", "-DNDEBUG"],
            extra_link_args=["-O3", "-fopenmp", "-DNDEBUG"]
        )
    ]
)
