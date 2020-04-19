from setuptools import setup, find_packages

setup(
    name='pytransport',
    version='1.3.1',
    packages=find_packages(exclude=["docs", "tests", "obsolete"]),
    # Not sure how strict these need to be...
    install_requires=["matplotlib",
                      "numpy",
                      "scipy"],
    # Some version of Python 3
    python_requires="==3.*",

    author='JAI@RHUL',
    author_email='william.shields.2010@live.rhul.ac.uk',
    description="Convert TRANSPORT models and load TRANSPORT output.",
    license='GPL3',
    url='https://bitbucket.org/jairhul/pytransport/'
)
