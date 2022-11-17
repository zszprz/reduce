# reduce

`reduce.py` is a python module that help with pairwise comparison matrices computations, such as CR reductions and calculate multiple indexes


### How to use it?

- It supports pairwise comparison matrices from 3x3 to 14x14
- Module `reduce.py` is in reduce folder of this repository
- Prepare CSV file as matrix.csv in main folder of this repository and use `import_pc_matrix_from_csv(filename)` or
- Generate random matrix by using `import_pc_matrix_from_csv(filename)`

### Consistency reduction algorithms:

<b>Xu and Wei:</b> `new_matrix, ci, cr = xu_and_wei_cr(matrix, lambda, threshold)`

- `matrix` is generated or imported matrix
- `lambda` is a lambda parameter for that algorithm (if you don't know what is it, set it to 0.9)
- `threshold` is a threshold to which value we want to reduce CR

<b>Cao:</b> `new_matrix, ci, cr = cao_cr(matrix, lambd, threshold)`

- `matrix` is generated or imported matrix
- `lambda` is a lambda parameter for that algorithm (if you don't know what is it, set it to 0.9)
- `threshold` is a threshold to which value we want to reduce CR

<b>Szybowski:</b> `new_matrix, ci, cr = szybowski_cr(matrix, threshold)`

- `matrix` is generated or imported matrix
- `threshold` is a threshold to which value we want to reduce CR

### Pairwise comparison matrices indexes:

<b>Koczkodaj Index (KI):</b> `index = koczkodaj_index(matrix)`

- `matrix` is generated or imported matrix

<b>Golden Wang Index (GWI):</b> `index = golden_wang_index(matrix)`

- `matrix` is generated or imported matrix

<b>Palaez Lamata Index (PLI):</b> `index = pelaez_lamata_index(matrix)`

- `matrix` is generated or imported matrix

<b>Geometric Consistency Index (GCI):</b> `index = geometric_consistency_index(matrix)`

- `matrix` is generated or imported matrix

<b>Triads Geometric Consistency Index (TGCI):</b> `index = triads_geometric_consistency_index(matrix)`

- `matrix` is generated or imported matrix

<b>Relative Error Index (REI):</b> `relative_error_index(matrix)`

- `matrix` is generated or imported matrix

<b>Harmonic Consistency Index (HCI):</b> `harmonic_consistency_index(matrix)`

- `matrix` is generated or imported matrix

### Development features

- Supports for `Python 3.9` and higher.
- [`Poetry`](https://python-poetry.org/) as the dependencies manager. See configuration in [`pyproject.toml`](https://github.com/zsz_prz/reduce/reduce/blob/master/pyproject.toml) and [`setup.cfg`](https://github.com/zsz_prz/reduce/reduce/blob/master/setup.cfg).
- Automatic codestyle with [`black`](https://github.com/psf/black), [`isort`](https://github.com/timothycrosley/isort) and [`pyupgrade`](https://github.com/asottile/pyupgrade).
- Ready-to-use [`pre-commit`](https://pre-commit.com/) hooks with code-formatting.
- Type checks with [`mypy`](https://mypy.readthedocs.io); docstring checks with [`darglint`](https://github.com/terrencepreilly/darglint); security checks with [`safety`](https://github.com/pyupio/safety) and [`bandit`](https://github.com/PyCQA/bandit)
- Testing with [`pytest`](https://docs.pytest.org/en/latest/).
- Ready-to-use [`.editorconfig`](https://github.com/zsz_prz/reduce/reduce/blob/master/.editorconfig), [`.dockerignore`](https://github.com/zsz_prz/reduce/reduce/blob/master/.dockerignore), and [`.gitignore`](https://github.com/zsz_prz/reduce/reduce/blob/master/.gitignore). You don't have to worry about those things.


## 🛡 License

This project is licensed under the terms of the `MIT` license. See [LICENSE](https://github.com/zsz_prz/reduce/reduce/blob/master/LICENSE) for more details.

## 📃 Citation

```bibtex
@misc{reduce,
  author = {ZSZ & SBA},
  title = {Reduce is a python package that help with pairwise comparison matrices computations, such as CR reductions and calculate multiple indexes},
  year = {2022},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/zsz_prz/reduce/}}
}
```

## Credits [![🚀 Your next Python package needs a bleeding-edge project structure.](https://img.shields.io/badge/python--package--template-%F0%9F%9A%80-brightgreen)](https://github.com/TezRomacH/python-package-template)

This project was generated with [`python-package-template`](https://github.com/TezRomacH/python-package-template)
