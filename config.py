"""config.py -- Python mirror of config.R.

Python scripts import the four directories from here:

    from config import RAW_DIR, DATA_DIR, RESULTS_DIR, FIG_DIR

Values can be overridden with the same environment variables that
config.R understands (RAW_DIR, DATA_DIR, RESULTS_DIR, FIG_DIR).
"""
import os

ROOT = os.path.dirname(os.path.abspath(__file__))

RAW_DIR     = os.environ.get("RAW_DIR",     os.path.join(ROOT, "data_raw"))
DATA_DIR    = os.environ.get("DATA_DIR",    os.path.join(ROOT, "data"))
RESULTS_DIR = os.environ.get("RESULTS_DIR", os.path.join(ROOT, "results"))
FIG_DIR     = os.environ.get("FIG_DIR",     os.path.join(ROOT, "figures"))

for _d in (DATA_DIR, RESULTS_DIR, FIG_DIR):
    os.makedirs(_d, exist_ok=True)

if __name__ == "__main__":
    print("ROOT       =", ROOT)
    print("RAW_DIR    =", RAW_DIR)
    print("DATA_DIR   =", DATA_DIR)
    print("RESULTS_DIR=", RESULTS_DIR)
    print("FIG_DIR    =", FIG_DIR)
