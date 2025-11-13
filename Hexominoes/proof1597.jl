include("Polyminoes.jl")
using .Polyminoes

board = parse(Board, read("reduced6.txt", String))
cp = canonical.(board.pieces)
freqs = Dict(x => count(==(x), cp) for x in unique(cp))
@info nholes(board)
omega = relax(cells(board); relax_by = 3)
optim(board, omega)

"""
This script demonstrates the key step in proving that 1597 is the
optimal hexomino configuration.

After performing the L-reduction with a penalty of 37, Gurobi is able to
establish local optimality in approximately 20 minutes.

The dots in proof1597.pdf represent the domain omega used in the
relaxation. Formally:

    There is no solution with more than 1634 = 1597 + 37 holes for this
    specific choice of L-shaped pieces that fits within the ring-shaped
    domain.

By enumerating and investigating all other possible L-shape reductions
(12 in total), this procedure establishes that 1597 is optimal for the max-fence
hexomino problem.
"""
