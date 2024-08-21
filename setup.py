import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SQDMetal",
    version="0.0.1.dev1",
    author="Prasanna Pakkiam",
    author_email="p.pakkiam@uq.edu.au",
    description="Tools to aid in simulating and fabricating superconducting quantum devices",
    url="https://github.com/sqdlab/SQDMetal",
    packages=setuptools.find_packages(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        "Programming Language :: Python :: 3 :: ONLY",
        "Operating System :: OS Independent",
    ],
    # python_requires='>=3.7',?
    keywords='superconducting qubits, comsol',
    install_requires=['qiskit-metal==0.1.2', 'mph', 'pyvista', 'jupyter', 'geopandas', 'gmsh'] #Qiskit-Metal doesn't install geopandas and Jupyter...
)
