from setuptools import setup, find_packages

setup(
    name = "RF-Tools",
    version = "0.0.0-dev",
    author = "John Garrett",
    author_email = "garrettj403@gmail.com",
    description = "RF tools",
    license = "MIT",
    url = "https://github.com/garrettj403/RF-tools/",
    packages=find_packages(),
    install_requires=[
        'matplotlib',
        'numpy',
        'scipy'
    ],
)
