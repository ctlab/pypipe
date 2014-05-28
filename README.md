pypipe
======

Bioinformatics pipelines framework,
writen in [python2](https://www.python.org/).

### Requirements
*   [python2](http://www.python.org/)
*   [setuptools](https://pythonhosted.org/setuptools/)
*   [graphviz](http://www.graphviz.org/)
*   [PyQt4](http://www.riverbankcomputing.co.uk/software/pyqt/intro) - If you want use GUI

### Instalation
    $ git clone https://github.com/semkagtn/pypipe.git
    $ cd pypipe
    $ python2 setup.py build
    $ sudo python2 setup.py install

### FAQ
*   **Installation for current user only:** `./setup.py install --user`
*   **How to write pipelines:** see example.py
*   **How to add support of new tool:** see pypipe/tools/\* files
*   **How to create pipeline:** `./pypipe.py create example.py pipeline`
*   **How to draw pipeline to img:** `./pypipe.py draw pipeline img.svg`
*   **How to run pipeline:** `./pypipe.py run pipeline 8`

