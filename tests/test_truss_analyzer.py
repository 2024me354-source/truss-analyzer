import csv
import os
import tempfile
import unittest
from unittest.mock import patch

from truss_analyzer import (
    TrussAnalyzer,
    compute_reactions,
    load_truss_from_file,
    load_truss_interactive,
    solve_joints,
    validate_truss,
)


WARREN = {
    "nodes": {"A": (0, 0), "B": (2, 0), "C": (4, 0), "D": (2, 2)},
    "members": [("A", "B"), ("B", "C"), ("A", "D"), ("B", "D"), ("C", "D")],
    "supports": {"A": "pin", "C": "roller"},
    "loads": {"D": (0, -10)},
}


class TrussAnalyzerTests(unittest.TestCase):
    def test_validate_truss_determinate(self):
        result = validate_truss(WARREN)
        self.assertTrue(result["is_determinate"])
        self.assertEqual(result["equation"], "5 + 3 = 8")

    def test_compute_reactions_warren(self):
        reactions = compute_reactions(WARREN)
        self.assertAlmostEqual(reactions["A"][0], 0.0, places=6)
        self.assertAlmostEqual(reactions["A"][1], 5.0, places=6)
        self.assertAlmostEqual(reactions["C"][0], 0.0, places=6)
        self.assertAlmostEqual(reactions["C"][1], 5.0, places=6)

    def test_solve_joints_warren(self):
        reactions = compute_reactions(WARREN)
        forces = solve_joints(WARREN, reactions)
        self.assertAlmostEqual(forces["AB"], 5.0, places=6)
        self.assertAlmostEqual(forces["BC"], 5.0, places=6)
        self.assertAlmostEqual(forces["AD"], -7.0710678119, places=6)
        self.assertAlmostEqual(forces["BD"], 0.0, places=6)
        self.assertAlmostEqual(forces["CD"], -7.0710678119, places=6)

    def test_truss_analyzer_end_to_end(self):
        analyzer = TrussAnalyzer(WARREN)
        forces = analyzer.solve()
        self.assertIn("AB", forces)
        table = analyzer.display_results()
        self.assertIn("TENSION", table)
        self.assertIn("COMPRESSION", table)

    def test_load_truss_from_csv(self):
        with tempfile.NamedTemporaryFile("w", newline="", suffix=".csv", delete=False) as f:
            writer = csv.writer(f)
            writer.writerow(["type", "node", "x", "y", "node1", "node2", "support", "fx", "fy"])
            writer.writerow(["node", "A", 0, 0, "", "", "", "", ""])
            writer.writerow(["node", "B", 2, 0, "", "", "", "", ""])
            writer.writerow(["member", "", "", "", "A", "B", "", "", ""])
            writer.writerow(["support", "A", "", "", "", "", "pin", "", ""])
            writer.writerow(["load", "B", "", "", "", "", "", 0, -5])
            f.flush()
        try:
            truss = load_truss_from_file(f.name)
            self.assertEqual(truss["nodes"]["A"], (0.0, 0.0))
            self.assertEqual(truss["members"][0], ("A", "B"))
            self.assertEqual(truss["supports"]["A"], "pin")
            self.assertEqual(truss["loads"]["B"], (0.0, -5.0))
        finally:
            os.unlink(f.name)

    def test_load_truss_from_csv_requires_type_column(self):
        with tempfile.NamedTemporaryFile("w", newline="", suffix=".csv", delete=False) as f:
            writer = csv.writer(f)
            writer.writerow(["node", "x", "y"])
            writer.writerow(["A", 0, 0])
            f.flush()
        try:
            with self.assertRaises(ValueError):
                load_truss_from_file(f.name)
        finally:
            os.unlink(f.name)

    def test_load_truss_interactive(self):
        answers = [
            "A", "0", "0",
            "B", "2", "0",
            "",  # end nodes
            "A", "B",
            "",  # end members
            "A", "pin",
            "",  # end supports
            "B", "0", "-5",
            "",  # end loads
        ]
        with patch("builtins.input", side_effect=answers):
            truss = load_truss_interactive()
        self.assertEqual(truss["nodes"]["A"], (0.0, 0.0))
        self.assertEqual(truss["members"][0], ("A", "B"))
        self.assertEqual(truss["supports"]["A"], "pin")
        self.assertEqual(truss["loads"]["B"], (0.0, -5.0))


if __name__ == "__main__":
    unittest.main()
