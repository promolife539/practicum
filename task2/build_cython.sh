#!/bin/bash

env CFLAGS="-I /Users/alex-gusev/.local/share/virtualenvs/task2-Lym8s-1u/lib/python3.6/site-packages/numpy/core/include $CFLAGS" pipenv run python setup.py build_ext --inplace