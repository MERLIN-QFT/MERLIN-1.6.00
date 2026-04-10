# MERLIN

**Method for Reduction of Loop Integrals**

MERLIN is a Mathematica-based program for the reduction of loop integrals in quantum field theory. It implements a covariant-derivative-based method to construct relations between Feynman integrals, providing an alternative to traditional reduction approaches.

---

## Abstract

We show how a large class of Feynman integrals can be efficiently reduced to master integrals by suitable covariant differentiation of the latter. The connections required for the covariant derivatives are constructed only once for a given topology and can then be applied to any configuration of internal propagator masses.

This algorithm is implemented in the Mathematica code MERLIN (Method for Reduction of Loop Integrals).

---

## Features

* Reduction of Feynman integrals to master integrals via covariant differentiation
* Applicability to arbitrary internal mass configurations
* Pre-built data for selected diagram families
* Symbolic implementation in Wolfram Mathematica

---

## Requirements

* Wolfram Mathematica (tested with version 14.3; version 12 or later recommended)

---

## Project Structure

MERLIN must be used with the provided directory structure. The working notebook (`.nb`) must be placed in the root directory, alongside the template notebook.

Typical structure:

```id="e8kq0r"
MERLIN/
├── data/
├── master-integrals/
├── matrices/
├── packages/
├── results/
└── MERLIN_TEMPLATE.nb
```

* `packages/` contains the core `.wl` implementation
* `data/`, `matrices/`, and `master-integrals/` store precomputed inputs
* `results/` is used for generated outputs
* The notebook must be executed from this root location to ensure correct path resolution

---

## Usage

MERLIN is used directly within a Mathematica notebook. No installation is required beyond maintaining the directory structure.

A minimal working example:

```mathematica id="9w3kzq"
SetDirectory[NotebookDirectory[]];
Get["packages/MERLIN.wl"];

INITIALIZE["1-loop-bubble"];
DIAGRAM[{2, 2}];
MASSCONFIG[u, u];

EVALUATE
```

### Workflow description

* `INITIALIZE[...]` initializes the chosen topology
* `DIAGRAM[...]` defines the diagram structure
* `MASSCONFIG[...]` sets the internal mass configuration
* `EVALUATE` performs the reduction

---

## Method

MERLIN implements a covariant derivative approach to generate relations among loop integrals. For a fixed topology, the required connections are constructed once and subsequently reused, allowing efficient reductions across different configurations of internal propagator masses.

---

## Pre-built Data

The current version includes precomputed data for the following diagram families:

* **1-loop bubble** → `"1-loop-bubble"`
* **1-loop triangle** → `"1-loop-triangle"`
* **2-loop vacuum bubble** → `"2-loops-vacuum"`
* **2-loop sunset** → `"2-loops-sunset"`
* **3-loop vacuum bubble** → `"3-loops-vacuum"`

These identifiers correspond to the arguments used in `INITIALIZE[...]`.

The data encode the necessary structures for the reduction procedure and can be extended by constructing additional topologies. (To be implemented)

Although some of these data were generated using external reduction tools, the construction is method-independent and can, in principle, be reproduced using other frameworks.

---

## Acknowledgments

Parts of the precomputed data included in this project were generated using FIRE (Feynman Integral REduction), which is licensed under the GNU GPL v2.0. Proper academic credit is given in the associated publication.

The construction of this data is not specific to FIRE and can, in principle, be reproduced using other reduction frameworks.

---

## License

This project is licensed under the GNU General Public License v2.0 (or later).

See the `LICENSE` file for details.

---

## Citation

If you use MERLIN in your research, please cite the associated publication:

> (Add reference here)

---

## Additional Notes

Users employing related tools (such as FIRE or LiteRed) should ensure proper citation of the corresponding publications, in accordance with their respective guidelines.

---

## Disclaimer

This software is provided without warranty of any kind and is intended for research use.
