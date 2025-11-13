include("Polyminoes.jl")
using .Polyminoes

pieces = generate(4)
omega = square_dom(7)
board = optim(pieces, omega; highs=true)
print(board)

"""

Finds

. . . . . . . .
. d d e e . . .
. d d . e e . .
. b . . . c . .
. b . . . c c .
. b b . . c . .
. . a a a a . .
. . . . . . . .

in a few seconds (with HiGHS instead of Gurobi)
"""


