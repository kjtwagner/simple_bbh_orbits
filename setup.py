from setuptools import find_packages, setup

setup(
    name='bbh',
    packages=find_packages(include=['bbh']),
    version='0.1.0',
    description='simple python library for bbh orbits',
    author='Kate Wagner',
    license='MIT',
    install_requires=['numpy','scipy']
)
