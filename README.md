# 2D Truss Analyzer

A comprehensive console-based 2D truss analyzer that solves trusses using the **Method of Joints**.

## Features

- **Static Determinacy Check**: Validates if truss is solvable (m + r = 2j)
- **Support Reaction Calculation**: Solves reactions using equilibrium equations (ΣFx=0, ΣFy=0, ΣM=0)
- **Method of Joints Solution**: Determines axial forces in all members
- **Force Classification**: Labels each member as TENSION, COMPRESSION, or ZERO FORCE
- **Results Table**: Clean formatted output of member forces
- **Visualization**: Matplotlib visualization with color-coded members (red=compression, blue=tension)
- **File Input**: Load truss definitions from .txt or .csv files

## Coming Soon

- Interactive user input
- File loading support
- Advanced visualization features

## Usage

```python
from truss_analyzer import TrussAnalyzer

# Define truss (Warren truss example)
truss_def = {
    'nodes': {
        'A': (0, 0),
        'B': (2, 0),
        'C': (4, 0),
        'D': (2, 2)
    },
    'members': [
        ('A', 'B'), ('B', 'C'), ('A', 'D'), ('B', 'D'), ('C', 'D')
    ],
    'supports': {
        'A': 'pin',
        'C': 'roller'
    },
    'loads': {
        'D': (0, -10)  # 10 kN downward
    }
}

# Solve truss
analyzer = TrussAnalyzer(truss_def)
analyzer.solve()
analyzer.display_results()
```

## Installation

```bash
pip install numpy matplotlib
```

## License

MIT
