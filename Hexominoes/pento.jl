include("Polyminoes.jl")
using .Polyminoes

pieces = generate(5)
omega = square_dom(19)
board = optim(pieces, omega)

"""
computes the optimal 128 solution for the pentomio problem

. . . . . . . . . . . . . . . . .
. . . . g . . . c . . . . . . . .
. . . . g g g c c c c f f f . . .
. . j j g . . . . . . . . f . . .
. j j . . . . . . . . . . f . . .
. j . . . . . . . . . . . h . . .
. a . . . . . . . . . . . h h . .
. a . . . . . . . . . . . . h . .
. a . . . . . . . . . . . . h . .
. a . . . . . . . . . . . . e e .
. a . . . . . . . . . . . . . e .
. b . . . . . . . . . . . . e e .
. b . . . . . . . . . . . k k . .
. b . . . . . . . . . . . k . . .
. b b i . . . l . . d d k k . . .
. . . i i i l l l d d d . . . . .
. . . . i . . l . . . . . . . . .
. . . . . . . . . . . . . . . . .

________________________________________________________
Executed in   28.11 mins    fish           external
   usr time   28.06 mins    1.13 millis   28.06 mins
   sys time    0.01 mins    0.00 millis    0.01 mins
"""

