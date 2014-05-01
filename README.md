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
*   **Installation for current user only:** python2 setup.py install --user
*   **How to write pipelines:** see example.py
*   **How to draw pipeline to png:** python2 pypipe.py example.py --draw
*   **How to run pipeline:** python2 pypipe.py example.py --run vcf
