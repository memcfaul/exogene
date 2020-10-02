#! /usr/bin/env python3#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Matt E. McFaul
@version: ExoGene 1.0.2
@email: memcfaul@ucdavis.edu
@copyright: Copyright 2017-2019, Draper Lab - UC Davis
@license: GPL 3.0+
"""

from xo_gui import exogene as xo
from tkinter import Tk
from pathlib import Path as plp
import sys

app_data = plp.home() / ".ExoGene"
if sys.platform == "darwin":
    app_data = plp.home() / "Library/Preferences/ExoGene"
elif sys.platform == "win32":
    app_data = plp.home() / "AppData/Roaming/ExoGene"
if not app_data.exists():
    app_data.mkdir()
with open(app_data / "log.txt", "w") as f:
    sys.stdout = f
    root = Tk()
    root.title("ExoGene 1.0.2")
    root.configure(background="#ececec")
    main = xo(root)
    root.protocol("WM_DELETE_WINDOW", main.quit)
    main.start()
    root.update()
    root.mainloop()
