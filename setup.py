import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="polyan-tobiasvonderhaar",
    version="2020.1",
    author="Tobias von der Haar",
    author_email="T.von-der-Haar@kent.ac.uk",
    description="A package for simulating polysome profiles from Ribo-Seq data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)