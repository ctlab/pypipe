pypipe
======

Bioinformatics pipelines framework,
Writen in python2

### Requirements
*   [python2](http://www.python.org/)
*   [setuptools](https://pythonhosted.org/setuptools/)
*   [graphviz](http://www.graphviz.org/)

### Instalation
    $ git clone https://github.com/semkagtn/pypipe.git
    $ cd pypipe
    $ python2 setup.py build
    $ sudo python2 setup.py install

### FAQ
*   **Installation for current user only:** ./setup.py install --user
*   **How to write pipelines:** see example.py
*   **How to add support of new tool:** see pypipe/tools/\* files
*   **How to draw pipeline to img:** ./pypipe.py example.py --draw svg
*   **How to run pipeline:** ./pypipe.py example.py --run vcf
