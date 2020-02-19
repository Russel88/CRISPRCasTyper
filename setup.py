import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="caspredict", 
    version="0.2.0",
    author="Jakob Russel",
    author_email="russel2620@gmail.com",
    description="Automatic detection and subtyping of CRISPR-Cas operons",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/russel88/caspredict",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    install_requires=[
        "numpy >= 1.17.5",
        "pandas >= 0.25.3",
        "scipy >= 1.4.1",
        "biopython >= 1.76",
        "multiprocess >= 0.70.9",
        "setuptools"],
    scripts=['bin/caspredict']
)
