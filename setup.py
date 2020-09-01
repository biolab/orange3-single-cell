#!/usr/bin/env python3
import sys
import setuptools
from os import path, walk


def include_documentation(local_dir, install_dir):
    data_files = []
    if (('bdist_wheel' in sys.argv) or ('install' in sys.argv))\
            and not path.exists(local_dir):
        print("Directory '{}' does not exist. "
              "Please build documentation before running bdist_wheel."
              .format(path.abspath(local_dir)))
        sys.exit(0)

    doc_files = []
    for dirpath, dirs, files in walk(local_dir):
        doc_files.append((dirpath.replace(local_dir, install_dir),
                          [path.join(dirpath, f) for f in files]))
    data_files.extend(doc_files)
    return data_files


if __name__ == '__main__':
    setuptools.setup(
        # data_files=include_documentation('doc/build/htmlhelp', 'help/orange3-single_cell'),
        use_scm_version=True,
        setup_requires=[
            'setuptools-scm',
            'setuptools',
        ],
        install_requires=[
            'Orange3>=3.23.0',
            'orange3-bioinformatics>=4.0.0',
            'fastdtw==0.3.2',
            'pandas>=0.23,<1.1',
            'loompy>=2.0.10',
            'xlrd~=1.2.0',
            'anndata>=0.6.21',
            'numpy',
            'scikit-learn',
        ],
        extras_require={
            'doc': ['sphinx', 'recommonmark'],
            'package': ['twine', 'wheel'],
            'test': [
                'pytest~=5.1.0',
                'pytest-cov~=2.7.1',
                'coverage~=4.5.4',
                'codecov~=2.0.15'
            ],
        },
    )
