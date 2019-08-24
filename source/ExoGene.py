from xo_gui import exogene as xo
from tkinter import Tk
from pathlib import Path as plp
import sys
output = plp.home() / "Library/Preferences/ExoGene/log.txt"
if not output.parents[0].exists():
    output.parents[0].mkdir()
with open(output, 'w') as f:
    sys.stdout = f
    root = Tk()
    root.geometry("555x357")
    root.update()
    root.title("ExoGene 1.0.1")
    main = xo(root)
    root.protocol("WM_DELETE_WINDOW", main.quit)
    main.start()
    root.update()
    root.geometry("558x359")
    root.mainloop()