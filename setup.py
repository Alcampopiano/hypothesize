from setuptools import find_packages, setup
import io
import os
import re

def get_install_requirements(path):
    content = read(path)
    return [req for req in content.split("\n") if req != "" and not req.startswith("#")]

def read(path, encoding="utf-8"):
    path = os.path.join(os.path.dirname(__file__), path)
    with io.open(path, encoding=encoding) as fp:
        return fp.read()

setup(
    name='hypothesize',
    version='0.1.dev13',
    description='A Python package for comparing groups and measuring associations using robust statistics.',
    author='Allan Campopiano',
    author_email="campopianoa@hcdsb.org",
    license='BSD 3-clause',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    url="https://github.com/Alcampopiano/hypothesize",
    packages=find_packages(),
    include_package_data=True,
    install_requires=get_install_requirements("requirements.txt"),
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 1 - Planning",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],

)
