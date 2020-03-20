import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="caspredict", 
    version="0.4.1",
    author="Jakob Russel",
    author_email="russel2620@gmail.com",
    description="Automatic detection and subtyping of CRISPR-Cas operons",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/russel88/caspredict",
    download_url="https://github.com/Russel88/CasPredict/archive/v0.4.1.tar.gz",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Development Status :: 4 - Beta"],
    python_requires='>=3.8',
    install_requires=[
        "numpy >= 1.17.5",
        "pandas >= 0.25.3",
        "scipy >= 1.4.1",
        "biopython >= 1.76",
        "multiprocess >= 0.70.9",
        "scikit-learn >= 0.22.0",
        "xgboost >= 0.9.0",
        "tqdm >= 4",
        "Pillow >= 7.0.0",
        "setuptools"],
    scripts=['bin/caspredict',
             'bin/repeatType']
)
