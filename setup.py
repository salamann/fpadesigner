import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="fpa",
    version="0.0.1",
    author="salamann",
    author_email="shunyo44+fpa@gmail.com",
    description="This application provides aerodynamics computation and flight dynamics of human powered airplane",
    license="MIT",
    keywords="aerodynamics flight human-powered",
    packages=["fpa"],
    long_description=read('README.md'),
    classifiers=[
    ],
)
