import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dnaseq-genbit",
    version="0.0.1",
    author="Akhilesh Balaji",
    description="Computationally operate on DNA and RNA sequences",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    py_modules=["genbit"],
    package_dir={'': 'src/genbit'},
    install_requires=["rich"]
)
