from setuptools import setup, find_packages

setup(
    name="flams",
    version="0.0.1",
    python_requires=">=3.10",
    packages=find_packages(
        include=["flams*"],
    ),
    install_requires=[
        "appdirs",
        "biopython",
        "requests",
    ],
)
