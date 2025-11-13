# The Pentomino Farm problem and other tesselations

## Overview

In ''The Colossal Book of Short Puzzles and Problems'' Martin Gardener gives a solution of 128 to the pentomino farm problem, where the 12 free pentominoes have to be placed edge-to-edge on the grid, so that they enclose the largest area possible. The solution was submitted by Donald Knuth and in 1978 Takakazu Shimauchi proved it was the maximum.

Here we present a Julia implementation for solving this problem using an ILP approach described by D.W. [here](https://cs.stackexchange.com/questions/153818/maximize-enclosed-area-of-given-figures-on-2d-grid/153851#153851). In addition to the original pentominoes, this repository extends its capabilities to examine additional tessellations. Specifically, it includes the analysis and solution for both tetrahexes and hexiamonds configurations. The implementation uses the commercial solver Gurobi and allows to either compute a single optimal solution or enumerate all possible solutions.

## Authors

### Original Pentomino Fence and Tessellations

* Alexis Langlois-Rémillard (alexislangloisremillard@gmail.com) https://alexisl-r.github.io/
* Mia Müßig (nienna@miamuessig.de) https://miamuessig.de/
* Erika Roldan (erika.roldan@ma.tum.de) https://www.erikaroldan.net/

### New Hexomino Fence

* Mykhailo Lyader
* Helmut Podhaisky (helmut.podhaisky@mathematik.uni-halle.de)

## License
This project is licensed under the MIT License - see LICENSE file for details. If you use this code for academic purposes, please cite the paper:

Langlois-Rémillard, Alexis, Müßig, Mia N. and Roldán, Érika. **Maximale Zäune mit Polyformen**, Mitteilungen der Deutschen Mathematiker-Vereinigung, vol. 33, no. 3, 2025, pp. 187-199. https://doi.org/10.1515/dmvm-2025-0056
