from setuptools import setup, find_packages

setup(
    name="hermessmiles",
    version="0.1.2",
    author="Srijit Seal",
    author_email="srijit@understanding.bio",
    description="SMILES comparison tool with overlap visualization",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/srijitseal/hermessmiles",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "rdkit",
        "Pillow",
        "molvs",
        "dimorphite-dl",
        "loguru",
        "pandarallel",
        "datetime",
        "pytest"
    ],
    extras_require={
        "test": ["pytest"]
    },
    entry_points={
        "console_scripts": [
            "hermessmiles=hermessmiles.cli:main"
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)