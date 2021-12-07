import setuptools

"""
Install locally:
>>> python setup.py install
"""

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


def read_requirements(fname):
    requirements = []
    with open(fname, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            requirements.append(line)


install_requires = ['xarray', 'pysolar', 'pyrttov', 'multiprocessing']

setuptools.setup(
    name="RTTOV-WRF",
    version="1.0.0",
    author="Lukas Kugler",
    author_email="lukas.kugler@univie.ac.at",
    description="Simulate SEVIRI satellite channels from WRF output",
    long_description=long_description,
    long_description_content_type="Markdown",
    url="https://github.com/lkugler/RTTOV-WRF",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7.3",
    install_requires=install_requires,
)
