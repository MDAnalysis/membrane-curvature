"""
MembraneCurvature MDAkit
A tool to calculate Membrane Curvature
"""
import sys
from setuptools import setup, find_packages
import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

# Check Python version
if sys.version_info[:2] < (3, 6):
    raise RuntimeError('Python version >=3.6 required.')

# Check numpy is installed
try:
    import numpy as np
except ImportError:
    print("MembraneCurvature requires NumPy.")


if __name__ == "__main__":

    with open("README.md", "r") as handle:
        long_description = handle.read()

    # set requirements
    install_requires = [
        'numpy>=1.16.0',
        'mdanalysis>=2.0.0b0',
        'mdanalysistests>=2.0.0',
    ]

setup(
    # Self-descriptive entries which should always be present
    name='membrane_curvature',
    author='Estefania Barreto-Ojeda',
    author_email='estefania.b.ojeda@gmail.com',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='MIT',

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=['numpy>=1.16.0'] + pytest_runner,

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # url='http://www.my_package.com',  # Website
    # install_requires=[],              # Required packages, pulls from pip if needed; do not use for Conda deployment
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    python_requires=">=3.6",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
