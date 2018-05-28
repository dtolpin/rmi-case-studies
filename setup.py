from setuptools import setup, find_packages
import rmi

setup(
    name="rmi",
    version=rmi.__version__,
    packages=find_packages(),

    # source code layout
    test_suite="rmi.test.test_suite",

    # Generating the command-line tool
    entry_points={
        "console_scripts": [
        ]
    },

    # dependencies
    install_requires=[],
    dependency_links=[]
)
