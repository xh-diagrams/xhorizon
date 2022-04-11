.. skeleton documentation master file, created by
   sphinx-quickstart on Fri Mar 25 12:56:13 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

skeleton
================================


Welcome!
--------

Hello and welcome to the skeleton package!

This package is a skeleton for Python packages allowing:
  - ``pip`` distribution and installation,
  - sphinx/readthedocs documentation.

See:
   - `GitHub <https://github.com/jcschindler01/skeleton>`_.
   - `Documetation <https://skeleton-jcschindler01.readthedocs.io/>`_.
   - `TestPyPi Distro <https://test.pypi.org/project/skeleton-JCSCHINDLER01/>`_.

Utilizes:
   - `Sphinx <https://www.sphinx-doc.org/en/master/>`_ documentation.

      - Local build docs with ``make html`` in docs folder. First ``make clean`` to regenerate toctrees.
      - ReadTheDocs should rebuild automatically on ``git push``.

   - `versioningit <https://versioningit.readthedocs.io/>`_ version single-sourcing in ``pyproject.toml``.

      - Versions are set by the git tags (or commit number if dirty).
      - On ``python -m build`` (or ``pip install -e .``) versioningit sets correct version in project metadata based on git data.
      - On ``import skeleton`` a ``__version__`` parameter is obtained using ``importlib.metadata``.

   - `APIdoc <https://www.sphinx-doc.org/en/master/man/sphinx-apidoc.html>`_ to auto-generate API reference.

      - Must manually run ``sphinx-apidoc`` to generate sphinx autodoc files. Automated by ``apidoc-rebuild.bat``.
      - Only need to regenerate when structure changes, not every update, because generates autodoc files.


To Develop:
   - Clone from Github.
   - Build and local install with ``pip install -e .`` in root folder.
   - Edit local source.
   - Access local version from anywhere with ``import skeleton`` in test scripts.
   - (Only need to rebuild if you want to update version metadata.)
   - Commit and push changes to remote as usual.

To Distribute:
   - Set version with ``git tag``, ``git push``, ``git push --tags``.
   - Clear ``dist`` folder of old versions.
   - Build with ``python -m build`` in root folder.
   - `Upload <https://packaging.python.org/en/latest/tutorials/packaging-projects/#uploading-the-distribution-archives>`_ with ``python -m twine upload --repository testpypi dist/*`` in root folder.
   - To test, ``pip install`` from appropriate repo inside a `venv <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment>`_.




Getting Started
---------------

If this were uploaded on PyPi (it's not, but you can download from TestPyPi link above), you would install the package with::

   pip install skeleton

Then, visit the tutorial to learn how to use it!


Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   Welcome! <welcome>
   Tutorial <TUTORIAL/00-overview>
   API Reference <API/apidoc/modules>


Search
------


* :ref:`Modules <modindex>`
* :ref:`Index <genindex>`
* :ref:`Search <search>`
