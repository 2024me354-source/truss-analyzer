# 2D Truss Analyzer

A comprehensive console-based 2D truss analyzer that solves trusses using the **Method of Joints**.

## Features

- Input from Python dictionary, interactive console prompts, or `.txt`/`.csv` files
- Static determinacy check (`m + r = 2j`)
- Support reaction calculation from global equilibrium (`ΣFx = 0`, `ΣFy = 0`, `ΣM = 0`)
- Joint-by-joint solution for member forces
- Member classification as **TENSION**, **COMPRESSION**, or **ZERO FORCE**
- Optional matplotlib drawing with blue tension and red compression members
- Graceful error handling for unstable/indeterminate trusses and solve failures

## Quick Start

Run the built-in Warren truss sample:

```bash
python truss_analyzer.py
```

Run with interactive input:

```bash
python truss_analyzer.py --interactive
```

Run from file input:

```bash
python truss_analyzer.py --file /path/to/truss.txt
python truss_analyzer.py --file /path/to/truss.csv
```

Optional visualization:

```bash
python truss_analyzer.py --plot
```

## Python API

```python
from truss_analyzer import TrussAnalyzer

truss_def = {
    "nodes": {
        "A": (0, 0),
        "B": (2, 0),
        "C": (4, 0),
        "D": (2, 2),
    },
    "members": [("A", "B"), ("B", "C"), ("A", "D"), ("B", "D"), ("C", "D")],
    "supports": {"A": "pin", "C": "roller"},
    "loads": {"D": (0, -10)},
}

analyzer = TrussAnalyzer(truss_def)
analyzer.solve()
analyzer.display_results()
```

## Installation

```bash
pip install numpy matplotlib
```

`numpy` and `matplotlib` are optional but recommended.

## License

MIT
