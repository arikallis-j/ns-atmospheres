from setuptools import setup, find_packages

# Creating packages for PyPi
#python3 setup.py bdist_wheel sdist
#twine upload dist/*

version = "0.1.0"
description = "Simple python library to work with neutron star's atmospheres"

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="atmons",
    version=version,
    description=description,
    package_dir={"atmons": "ns","atmons": "cli"},
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/arikallis-j/ns-atmospheres",
    author="Arik Allis",
    author_email="arikallis.j@gmail.com",
    license="GPL-2.0",
    classifiers=[
    "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    "Programming Language :: Python :: 3.10",
    "Operating System :: OS Independent",
    "Development Status :: 3 - Alpha"
    "Topic :: Scientific/Engineering :: Astronomy",
    "Intended Audience :: Science/Research",
    "Natural Language :: English"
    ],
    packages=find_packages(),
    install_requires=[
        "yaml",
        "numpy",
        "astropy",
        "matplotlib",
        "rich",
    ],
    python_requires=">=3.8",
)