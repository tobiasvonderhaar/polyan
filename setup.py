import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="polyan",
    version="2020.1.2",
    author="Tobias von der Haar",
    author_email="T.von-der-Haar@kent.ac.uk",
    license = "MIT",
    description="A package for simulating polysome profiles from Ribo-Seq data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tobiasvonderhaar/polyan",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
        "Framework :: Jupyter",
        "Framework :: Matplotlib",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    python_requires='>=3.6',
    install_requires=[
          'pandas','numpy','matplotlib','scipy'
      ],
)