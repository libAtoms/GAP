import setuptools

setuptools.setup(
    name="hybrid_md",
    version="0.0.1",
    packages=setuptools.find_packages(),
    install_requires=["click>=7.0", "numpy", "ase", "pyyaml"],
    entry_points="""
    [console_scripts]
    wfl=hybrid_md.cli:main
    """,
)
