# https://packaging.python.org/tutorials/packaging-projects/
import versioneer
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="objectgmsh",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Arved Enders-Seidlitz",
    author_email="arved.enders-seidlitz@ikz-berlin.de",
    description="Object oriented Gmsh modeling.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nemocrys/objectgmsh",
    packages=["objectgmsh", "objectgmsh.test"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "gmsh",
    ],
    python_requires=">=3.7",
)
