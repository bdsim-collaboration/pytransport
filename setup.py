from setuptools import setup, find_packages

setup(
    name='pytransport',
    version='1.5.0',
    packages=find_packages(exclude=["docs", "tests"]),
    # Not sure how strict these need to be...
    install_requires=["matplotlib",
                      "numpy",
                      "root_numpy",
                      "scipy"],
    python_requires=">=3.*",

    author='JAI@RHUL',
    author_email='william.shields@rhul.ac.uk',
    description="Convert TRANSPORT models and load TRANSPORT output.",
    license='GPL3',
    url='https://bitbucket.org/jairhul/pytransport/'
)
