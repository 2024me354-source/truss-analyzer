"""Console-based 2D truss analyzer using the Method of Joints."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

try:
    import numpy as np  # type: ignore
except ImportError:  # pragma: no cover - optional dependency
    np = None


Vector = Tuple[float, float]
Member = Tuple[str, str]
SOLVE_TOLERANCE = 1e-9
EQUILIBRIUM_TOLERANCE = 1e-6
GEOMETRIC_TOLERANCE = 1e-12


@dataclass
class TrussData:
    nodes: Dict[str, Vector]
    members: List[Member]
    supports: Dict[str, str]
    loads: Dict[str, Vector]


def _normalize_member(member: Iterable[str]) -> Member:
    a, b = tuple(member)
    if a == b:
        raise ValueError(f"Member endpoints must differ: {a}-{b}")
    return str(a), str(b)


def _member_name(a: str, b: str) -> str:
    if len(a) == 1 and len(b) == 1:
        return f"{a}{b}"
    return f"{a}-{b}"


def _as_truss_data(truss_def: Dict) -> TrussData:
    if not isinstance(truss_def, dict):
        raise TypeError("Truss definition must be a dictionary.")

    nodes_raw = truss_def.get("nodes", {})
    members_raw = truss_def.get("members", [])
    supports_raw = truss_def.get("supports", {})
    loads_raw = truss_def.get("loads", {})

    nodes = {str(k): (float(v[0]), float(v[1])) for k, v in nodes_raw.items()}
    members = [_normalize_member(m) for m in members_raw]
    supports = {str(k): str(v).lower() for k, v in supports_raw.items()}
    loads = {str(k): (float(v[0]), float(v[1])) for k, v in loads_raw.items()}

    if not nodes:
        raise ValueError("No nodes defined.")
    if not members:
        raise ValueError("No members defined.")

    for a, b in members:
        if a not in nodes or b not in nodes:
            raise ValueError(f"Member {_member_name(a, b)} references undefined node.")

    for node, support_type in supports.items():
        if node not in nodes:
            raise ValueError(f"Support node '{node}' is undefined.")
        if support_type not in {"pin", "roller"}:
            raise ValueError(f"Unsupported support type '{support_type}' at {node}.")

    for node in loads:
        if node not in nodes:
            raise ValueError(f"Load node '{node}' is undefined.")

    return TrussData(nodes=nodes, members=members, supports=supports, loads=loads)


def validate_truss(truss_def: Dict) -> Dict[str, object]:
    """Validate static determinacy using m + r = 2j."""
    truss = _as_truss_data(truss_def)
    m = len(truss.members)
    r = 0
    for support_type in truss.supports.values():
        r += 2 if support_type == "pin" else 1
    j = len(truss.nodes)

    lhs = m + r
    rhs = 2 * j
    if lhs == rhs:
        status = "statically determinate"
    elif lhs > rhs:
        status = "statically indeterminate"
    else:
        status = "mechanism (unstable)"

    return {
        "m": m,
        "r": r,
        "j": j,
        "equation": f"{m} + {r} = {2*j}",
        "status": status,
        "is_determinate": lhs == rhs,
    }


def compute_reactions(truss_def: Dict) -> Dict[str, Vector]:
    """Compute support reactions from global equilibrium equations."""
    truss = _as_truss_data(truss_def)

    pin_nodes = [n for n, t in truss.supports.items() if t == "pin"]
    roller_nodes = [n for n, t in truss.supports.items() if t == "roller"]

    if len(pin_nodes) != 1 or len(roller_nodes) != 1:
        raise ValueError("Current solver supports exactly one pin and one roller.")

    pin = pin_nodes[0]
    roller = roller_nodes[0]

    x_pin, y_pin = truss.nodes[pin]
    x_roller, _ = truss.nodes[roller]
    lever = x_roller - x_pin
    if abs(lever) < GEOMETRIC_TOLERANCE:
        raise ValueError("Pin and roller cannot be vertically aligned for this solver.")

    sum_fx = sum(force[0] for force in truss.loads.values())
    sum_fy = sum(force[1] for force in truss.loads.values())

    moment_loads_about_pin = 0.0
    for node, (fx, fy) in truss.loads.items():
        x, y = truss.nodes[node]
        moment_loads_about_pin += (x - x_pin) * fy - (y - y_pin) * fx

    roller_ry = -moment_loads_about_pin / lever
    pin_ry = -sum_fy - roller_ry
    pin_rx = -sum_fx

    reactions = {node: (0.0, 0.0) for node in truss.nodes}
    reactions[pin] = (pin_rx, pin_ry)
    reactions[roller] = (0.0, roller_ry)
    return reactions


def _unit_vector(a: Vector, b: Vector) -> Vector:
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    length = math.hypot(dx, dy)
    if length == 0:
        raise ValueError("Zero-length member detected.")
    return dx / length, dy / length


def solve_joints(truss_def: Dict, reactions: Optional[Dict[str, Vector]] = None) -> Dict[str, float]:
    """Solve member axial forces by Method of Joints.

    Positive value means tension, negative means compression.
    """
    truss = _as_truss_data(truss_def)
    reactions = reactions or compute_reactions(truss_def)

    member_forces: Dict[Member, Optional[float]] = {m: None for m in truss.members}
    connected: Dict[str, List[Member]] = {node: [] for node in truss.nodes}
    for m in truss.members:
        connected[m[0]].append(m)
        connected[m[1]].append(m)

    progress = True
    while progress:
        progress = False
        for joint in truss.nodes:
            unknowns: List[Tuple[Member, float, float]] = []
            known_x = truss.loads.get(joint, (0.0, 0.0))[0] + reactions.get(joint, (0.0, 0.0))[0]
            known_y = truss.loads.get(joint, (0.0, 0.0))[1] + reactions.get(joint, (0.0, 0.0))[1]

            for m in connected[joint]:
                a, b = m
                ux, uy = _unit_vector(truss.nodes[a], truss.nodes[b])
                cx, cy = (ux, uy) if joint == a else (-ux, -uy)

                force = member_forces[m]
                if force is None:
                    unknowns.append((m, cx, cy))
                else:
                    known_x += cx * force
                    known_y += cy * force

            if not unknowns or len(unknowns) > 2:
                continue

            rhs_x = -known_x
            rhs_y = -known_y

            if len(unknowns) == 1:
                member, c1x, c1y = unknowns[0]
                if abs(c1x) > SOLVE_TOLERANCE:
                    force_value = rhs_x / c1x
                elif abs(c1y) > SOLVE_TOLERANCE:
                    force_value = rhs_y / c1y
                else:
                    continue
                if (
                    abs(c1x * force_value - rhs_x) > EQUILIBRIUM_TOLERANCE
                    or abs(c1y * force_value - rhs_y) > EQUILIBRIUM_TOLERANCE
                ):
                    continue
                member_forces[member] = force_value
                progress = True
                continue

            (m1, c1x, c1y), (m2, c2x, c2y) = unknowns
            det = c1x * c2y - c1y * c2x
            if abs(det) < SOLVE_TOLERANCE:
                if np is None:
                    continue
                matrix = np.array([[c1x, c2x], [c1y, c2y]], dtype=float)
                vec = np.array([rhs_x, rhs_y], dtype=float)
                try:
                    f1, f2 = np.linalg.solve(matrix, vec).tolist()
                except np.linalg.LinAlgError:
                    continue
            else:
                f1 = (rhs_x * c2y - rhs_y * c2x) / det
                f2 = (c1x * rhs_y - c1y * rhs_x) / det

            member_forces[m1] = f1
            member_forces[m2] = f2
            progress = True

    unresolved = [_member_name(a, b) for (a, b), f in member_forces.items() if f is None]
    if unresolved:
        raise ValueError(
            f"Unable to solve all member forces with Method of Joints. Unresolved members: "
            f"{', '.join(unresolved)}. This may indicate an unstable configuration, improper "
            f"support constraints, or members not connected in a solvable sequence."
        )

    return {_member_name(a, b): float(member_forces[(a, b)]) for a, b in truss.members}


def classify_force(force: float, tolerance: float = SOLVE_TOLERANCE) -> str:
    if abs(force) <= tolerance:
        return "ZERO FORCE"
    return "TENSION" if force > 0 else "COMPRESSION"


def display_results(member_forces: Dict[str, float], unit: str = "kN") -> str:
    """Return and print a formatted member force table."""
    header = ["Member", f"Force ({unit})", "Nature"]
    line = "+------------+--------------+-------------+"
    rows = [line, f"| {header[0]:<10} | {header[1]:<12} | {header[2]:<11} |", line]

    for member in sorted(member_forces):
        force = member_forces[member]
        nature = classify_force(force)
        rows.append(f"| {member:<10} | {abs(force):>12.4f} | {nature:<11} |")

    rows.append(line)
    table = "\n".join(rows)
    print(table)
    return table


def draw_truss(
    truss_def: Dict,
    member_forces: Optional[Dict[str, float]] = None,
    title: str = "2D Truss Analyzer",
) -> None:
    """Draw truss with optional force coloring (blue=tension, red=compression)."""
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("matplotlib is required for visualization") from exc

    truss = _as_truss_data(truss_def)
    fig, ax = plt.subplots()

    for a, b in truss.members:
        xa, ya = truss.nodes[a]
        xb, yb = truss.nodes[b]
        name = _member_name(a, b)
        color = "black"
        if member_forces is not None and name in member_forces:
            nature = classify_force(member_forces[name], tolerance=SOLVE_TOLERANCE)
            color = "blue" if nature == "TENSION" else "red" if nature == "COMPRESSION" else "gray"
        ax.plot([xa, xb], [ya, yb], color=color, linewidth=2)
        ax.text((xa + xb) / 2, (ya + yb) / 2, name, fontsize=9)

    for node, (x, y) in truss.nodes.items():
        ax.plot(x, y, "ko")
        ax.text(x + 0.05, y + 0.05, node)

    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title)
    ax.grid(True, alpha=0.2)
    plt.show()


def load_truss_from_file(path: str) -> Dict:
    """Load truss from .txt or .csv file."""
    if path.lower().endswith(".csv"):
        return _load_from_csv(path)
    if path.lower().endswith(".txt"):
        return _load_from_txt(path)
    raise ValueError("Only .txt and .csv files are supported.")


def _load_from_csv(path: str) -> Dict:
    nodes: Dict[str, Vector] = {}
    members: List[Member] = []
    supports: Dict[str, str] = {}
    loads: Dict[str, Vector] = {}

    with open(path, newline="", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        if not reader.fieldnames or "type" not in {f.lower() for f in reader.fieldnames}:
            raise ValueError(
                "CSV must include a 'type' column with values: node, member, support, or load."
            )

        field_map = {f.lower(): f for f in reader.fieldnames}
        for row in reader:
            row_type = row.get(field_map["type"], "").strip().lower()
            if row_type == "node":
                node = row.get(field_map.get("node", "node"), "").strip()
                nodes[node] = (
                    float(row.get(field_map.get("x", "x"), 0.0)),
                    float(row.get(field_map.get("y", "y"), 0.0)),
                )
            elif row_type == "member":
                n1 = row.get(field_map.get("node1", "node1"), "").strip()
                n2 = row.get(field_map.get("node2", "node2"), "").strip()
                members.append((n1, n2))
            elif row_type == "support":
                node = row.get(field_map.get("node", "node"), "").strip()
                supports[node] = row.get(field_map.get("support", "support"), "").strip().lower()
            elif row_type == "load":
                node = row.get(field_map.get("node", "node"), "").strip()
                loads[node] = (
                    float(row.get(field_map.get("fx", "fx"), 0.0)),
                    float(row.get(field_map.get("fy", "fy"), 0.0)),
                )

    return {"nodes": nodes, "members": members, "supports": supports, "loads": loads}


def _load_from_txt(path: str) -> Dict:
    nodes: Dict[str, Vector] = {}
    members: List[Member] = []
    supports: Dict[str, str] = {}
    loads: Dict[str, Vector] = {}

    section = ""
    with open(path, encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            normalized = line.strip("[]:").upper()
            if normalized in {"NODES", "MEMBERS", "SUPPORTS", "LOADS"}:
                section = normalized
                continue

            parts = [p.strip() for p in line.split(",")]
            if section == "NODES" and len(parts) == 3:
                nodes[parts[0]] = (float(parts[1]), float(parts[2]))
            elif section == "MEMBERS" and len(parts) == 2:
                members.append((parts[0], parts[1]))
            elif section == "SUPPORTS" and len(parts) == 2:
                supports[parts[0]] = parts[1].lower()
            elif section == "LOADS" and len(parts) == 3:
                loads[parts[0]] = (float(parts[1]), float(parts[2]))

    return {"nodes": nodes, "members": members, "supports": supports, "loads": loads}


def load_truss_interactive() -> Dict:
    """Collect truss data via console prompts."""
    print("Enter node data (name x y), blank name to stop.")
    nodes: Dict[str, Vector] = {}
    while True:
        name = input("Node name: ").strip()
        if not name:
            break
        x = float(input(f"{name} x: ").strip())
        y = float(input(f"{name} y: ").strip())
        nodes[name] = (x, y)

    print("Enter members (node1 node2), blank node1 to stop.")
    members: List[Member] = []
    while True:
        n1 = input("Member node1: ").strip()
        if not n1:
            break
        n2 = input("Member node2: ").strip()
        members.append((n1, n2))

    print("Enter supports (node type[pin/roller]), blank node to stop.")
    supports: Dict[str, str] = {}
    while True:
        node = input("Support node: ").strip()
        if not node:
            break
        supports[node] = input("Support type: ").strip().lower()

    print("Enter loads (node fx fy), blank node to stop.")
    loads: Dict[str, Vector] = {}
    while True:
        node = input("Load node: ").strip()
        if not node:
            break
        fx = float(input("Fx: ").strip())
        fy = float(input("Fy: ").strip())
        loads[node] = (fx, fy)

    return {"nodes": nodes, "members": members, "supports": supports, "loads": loads}


class TrussAnalyzer:
    def __init__(self, truss_def: Dict):
        self.truss_def = truss_def
        self.validation: Optional[Dict[str, object]] = None
        self.reactions: Optional[Dict[str, Vector]] = None
        self.member_forces: Optional[Dict[str, float]] = None

    def solve(self) -> Dict[str, float]:
        self.validation = validate_truss(self.truss_def)
        if not self.validation["is_determinate"]:
            raise ValueError(f"Truss is {self.validation['status']}: {self.validation['equation']}")
        self.reactions = compute_reactions(self.truss_def)
        self.member_forces = solve_joints(self.truss_def, self.reactions)
        return self.member_forces

    def display_results(self, unit: str = "kN") -> str:
        if self.member_forces is None:
            raise ValueError("Call solve() before display_results().")
        return display_results(self.member_forces, unit=unit)

    def draw_truss(self) -> None:
        draw_truss(self.truss_def, self.member_forces)


def _default_warren_truss() -> Dict:
    return {
        "nodes": {"A": (0, 0), "B": (2, 0), "C": (4, 0), "D": (2, 2)},
        "members": [("A", "B"), ("B", "C"), ("A", "D"), ("B", "D"), ("C", "D")],
        "supports": {"A": "pin", "C": "roller"},
        "loads": {"D": (0, -10)},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Console-based 2D Truss Analyzer")
    parser.add_argument("--file", help="Load truss from .txt or .csv")
    parser.add_argument("--interactive", action="store_true", help="Enter truss data interactively")
    parser.add_argument("--plot", action="store_true", help="Plot truss using matplotlib")
    parser.add_argument("--unit", default="kN", help="Display unit label")
    args = parser.parse_args()

    if args.file:
        truss_def = load_truss_from_file(args.file)
    elif args.interactive:
        truss_def = load_truss_interactive()
    else:
        truss_def = _default_warren_truss()

    analyzer = TrussAnalyzer(truss_def)
    validation = validate_truss(truss_def)
    print(
        f"Determinacy check: m+r=2j -> {validation['equation']} ({validation['status']})"
    )
    if not validation["is_determinate"]:
        return

    analyzer.solve()
    if analyzer.reactions:
        print("\nSupport Reactions:")
        for node, (rx, ry) in analyzer.reactions.items():
            if abs(rx) > GEOMETRIC_TOLERANCE or abs(ry) > GEOMETRIC_TOLERANCE:
                print(f"  {node}: Rx={rx:.4f}, Ry={ry:.4f}")
    print()
    analyzer.display_results(unit=args.unit)

    if args.plot:
        analyzer.draw_truss()


if __name__ == "__main__":
    main()
