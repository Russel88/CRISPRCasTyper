import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cctyper", 
    version="1.5.0",
    author="Jakob Russel",
    author_email="russel2620@gmail.com",
    description="CRISPRCasTyper: Automatic detection and subtyping of CRISPR-Cas operons",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Russel88/CRISPRCasTyper",
    download_url="https://github.com/Russel88/CRISPRCasTyper/archive/v1.5.0.tar.gz",
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
        "pandas >= 1.3",
        "scipy >= 1.4.1",
        "biopython >= 1.76",
        "multiprocess >= 0.70.9",
        "scikit-learn >= 0.22.0",
        "xgboost >= 1.4",
        "tqdm >= 4",
        "drawSvg >= 1.8.0",
        "setuptools"],
    scripts=['bin/cctyper',
             'bin/repeatType',
             'bin/repeatTrain']
)
