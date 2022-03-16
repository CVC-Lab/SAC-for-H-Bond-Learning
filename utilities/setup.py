from distutils.core import setup

setup(
    name='utilities',
    version='0.1',
    description='utility library',
    author='Conrad Li',
    author_email='conradliste@utexas.edu',
    requires=[ 'numpy', 'matplotlib',],
    packages=['utilities', ],
    package_data={
    	'utilities': ['*'],
    	'utilities.utils': ['*'],
    },
)