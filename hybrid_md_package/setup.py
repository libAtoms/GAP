#  Hybrid MD decision making package
#
#  Copyright (c) Tamas K. Stenczel 2021.
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hybrid_md",
    version="0.0.1",
    packages=setuptools.find_packages(),
    install_requires=["click>=7.0", "numpy", "ase", "pyyaml"],
    entry_points="""
    [console_scripts]
    hybrid-md=hybrid_md.cli:main
    """,
    description="Hybrid Molecular Dynamics for Quantum Mechanics codes with Force Fields",
    author="Tamas K. Stenczel",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
    ],
    zip_safe=True,
)
