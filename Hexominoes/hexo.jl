include("Polyminoes.jl")
using .Polyminoes

sx = 47
k = 6
m = sx/2
omega = Set((x, y) for (x, y) in circle_dom(sx) if (x-m)^2+(y-m)^2>(m-3)^2)

pieces = generate(k)
board = optim(pieces, omega)
println(board)
board = stochopt(board)
println(board)

"""
Computes a feasible initial configuration with 1496 holes and then
optimises that randomly with local relaxation.
In about 15 minutes a configuration with typically between
1590 to 1597 holes is obained (in about 15 minutes).
"""
