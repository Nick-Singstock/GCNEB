#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: zaba1157
"""

import os
from pathlib import Path

scripts = [script for script in os.listdir('.') if os.path.isfile(script) and Path(script).suffix == '.py']

for s in scripts:
    os.system('chmod 755 '+s)

# makes all python scripts (including setup.py) in the present working directory executable
# run: python setup.py
