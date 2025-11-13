"""
This code solves the max fence problem for various polyminoes configuration.
It is based on <https://github.com/PhoenixSmaug/Pentomino-Farm/blob/main/pentomino.jl>. 
It also uses local adaptivity, a stochastic heuristic and a reduction to L shapes.

See `tetro.jl`, `pento.jl`, `hexo.jl` and `proof1597.jl` how it can be used. 

Copyright (c) 2025 Mykhailo Lyader,
              Helmut Podhaisky <helmut.podhaisky@mathematik.uni-halle.de>    
"""
module Polyminoes
using JuMP, Gurobi, HiGHS

export Polymino,
    Pt,
    CellSet,
    Board,
    generate,
    hull,
    bbox,
    canonical,
    normalize,
    optim,
    cells,
    pieces,
    relax,
    holes,
    nholes,
    square_dom,
    circle_dom,
    n4,
    n8,
    findchain,
    stochopt

const Pt = Tuple{Int,Int}
const CellSet = Set{Pt}

function square_dom(sx::Int)::CellSet
    Set((x, y) for x = 1:sx, y = 1:sx)
end

function circle_dom(sx::Int)::CellSet
    m = sx/2
    Set((x, y) for (x, y) in square_dom(sx) if (x-m)^2 + (y-m)^2 <= m^2)
end


struct Polymino
    cells::Vector{Pt}
    function Polymino(cells::Vector{Pt})
        new(sort(cells))
    end
end

struct Board
    pieces::Vector{Polymino}
    function Board(pieces::Vector{Polymino})
        new(sort(pieces))
    end
end

Base.iterate(b::Board, state...) = iterate(b.pieces, state...)
Base.length(b::Board) = length(b.pieces)
Base.eltype(::Type{Board}) = Polymino
Base.getindex(b::Board, i::Int) = b.pieces[i]


function Base.show(io::IO, board::Board)
    ((x0, x1), (y0, y1)) = bbox(board)
    m = fill(". ", y1 - y0 + 3, x1 - x0 + 3)

    letters = vcat(
        [string(c, ' ') for c = 'a':'z'],
        [string(c, ' ') for c = 'A':'Z'],
        [string(c, d) for c = 'a':'z' for d = 'a':'z'],
    )

    bpieces = sort(board.pieces, by = canonical)

    for (lb, p) in zip(letters, bpieces)
        for xy in p.cells
            m[y1-xy[2]+2, 2+xy[1]-x0] = lb
        end
    end

    print(io, join([join(row, "") for row in eachrow(m)], "\n"))
end

function Base.parse(::Type{Board}, ascii::AbstractString)
    lines = split(ascii, '\n', keepempty = false)
    height = length(lines)
    width = maximum(length.(lines))

    grid = [
        begin
            lj = lines[j]
            if isodd(length(lj))
                lj = lj * " "
            end
            i <= length(lj) ? lj[i:(i+1)] : "  "
        end for j = 1:height, i = 1:2:width
    ]

    letters = Set(filter(s -> length(s) == 2 && isletter(s[1]), vec(grid)))

    pieces = Polymino[]
    for code in letters
        cells = Pt[]
        for j = 1:height, i = 1:(width÷2)
            if grid[j, i] == code
                push!(cells, (i, height - j + 1))
            end
        end
        p = Polymino(cells)
        push!(pieces, p)
    end

    Board(pieces)
end



Base.show(io::IO, p::Polymino) = print(io, "Polymino(", p.cells, ")")

const D8 = [
    (x, y) -> (x, y),
    (x, y) -> (-x, y),
    (x, y) -> (x, -y),
    (x, y) -> (-x, -y),
    (x, y) -> (y, x),
    (x, y) -> (-y, x),
    (x, y) -> (y, -x),
    (x, y) -> (-y, -x),
]

normalize(cells) = begin
    xs = map(first, cells)
    ys = map(last, cells)
    dx, dy = minimum(xs), minimum(ys)
    sort([(x - dx, y - dy) for (x, y) in cells])
end

function normalize(p::Polymino)::Polymino
    Polymino(normalize(p.cells))
end

function canonical(p::Polymino)::Polymino
    Polymino(minimum(normalize([f(x, y) for (x, y) in p.cells]) for f in D8))
end

Base.hash(p::Polymino, h::UInt) = hash(p.cells, h)
Base.:(==)(a::Polymino, b::Polymino) = a.cells == b.cells
Base.isless(a::Polymino, b::Polymino) = a.cells < b.cells

"""
Orbit of p under the dieder group.
"""
function sym(p::Polymino)
    variants = Set{Vector{Tuple{Int,Int}}}()
    for f in D8
        transformed = [f(x, y) for (x, y) in p.cells]
        push!(variants, normalize(transformed))
    end
    sorted_variants = sort(collect(variants))
    [Polymino(v) for v in sorted_variants]
end

function cells(p::Polymino)::CellSet
    CellSet(p.cells)
end

"""
Extend p by one piece
"""
function extend_polyomino(p::Polymino, k::Int, seen::Set{Polymino})
    if length(p.cells) == k
        push!(seen, canonical(p))
        return
    end
    for xy in p.cells
        for nb in n4(xy)
            if !(nb in p.cells)
                newp = Polymino(vcat(p.cells, [nb]))
                extend_polyomino(newp, k, seen)
            end
        end
    end
end

"""
All polyminos of size k
"""
function generate(k::Int)
    seen = Set{Polymino}()
    extend_polyomino(Polymino([(0, 0)]), k, seen)
    sort(collect(seen))
end


"""
Lift n4 to a polynomino p
"""
function hull(p::Polymino)::CellSet
    setdiff(CellSet([z for l in cells(p) for z in n4(l)]), cells(p))
end

"""
All ways to put piece p into domain omega.
"""
function place_at(xy::Pt, p::Polymino, omega::CellSet, rotate::Bool)::Vector{Polymino}
    variants = rotate ? sym(p) : [p]
    res = Polymino[]
    for q in variants
        for (dx, dy) in q.cells
            shift = (xy[1] - dx, xy[2] - dy)
            placed = Polymino([(x + shift[1], y + shift[2]) for (x, y) in q.cells])
            if all(c -> c in omega, placed.cells)
                push!(res, placed)
            end
        end
    end
    res
end

"""
The number of holes
"""
function nholes(board::Board)::Int
    length(holes(board))
end

function pieces(board::Board)::Vector{Polymino}
    collect(board.pieces)
end

function holes(board::Board)::CellSet
    holes(cells(board))
end

function holes(occupied::CellSet)::CellSet

    ((x0, x1), (y0, y1)) = bbox(occupied)

    filled = CellSet()

    queue = union(
        [(x, y) for x = x0:x1, y in (y0-1, y1+1)],
        [(x, y) for x in (x0-1, x1+1), y = y0:y1],
    )
    queue = collect(filter(c -> !(c in occupied), queue))
    while !isempty(queue)
        c = pop!(queue)
        if c in filled
            continue
        end
        push!(filled, c)
        for (nx, ny) in n4(c)
            if x0 <= nx <= x1 && y0 <= ny <= y1 && (nx, ny) ∉ occupied && (nx, ny) ∉ filled
                push!(queue, (nx, ny))
            end
        end
    end
    Set((x, y) for x in x0:x1, y in y0:y1 if !((x, y) in occupied) && !((x, y) in filled))
end

function cells(board::Board)::CellSet
    reduce(vcat, getfield.(board.pieces, :cells)) |> CellSet
end

function bbox(board::Board)::Tuple{Tuple{Int,Int},Tuple{Int,Int}}
    bbox(cells(board))
end

function bbox(all_cells::CellSet)::Tuple{Tuple{Int,Int},Tuple{Int,Int}}
    min_x, max_x = minimum(getindex.(all_cells, 1)), maximum(getindex.(all_cells, 1))
    min_y, max_y = minimum(getindex.(all_cells, 2)), maximum(getindex.(all_cells, 2))
    ((min_x, max_x), (min_y, max_y))
end

function grow(omega::CellSet)::CellSet
    CellSet(y for x in omega for y in n8(x))
end

function relax(omega::CellSet; relax_by = 2)::CellSet
    foldl((w, _) -> grow(w), 1:relax_by; init = omega)
end


function n8(xy::Pt)::Vector{Pt}
    x, y = xy
    [(x+dx, y+dy) for dx in (-1, 0, 1), dy in (-1, 0, 1) if (dx, dy) != (0, 0)]
end

function n4(xy::Pt)::Vector{Pt}
    x, y = xy
    [(x+dx, y+dy) for (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1))]
end

"""
Solves the max fence problem with the pieces.
Pieces outside omega are fixed in their position.
The indicator constraint connect2 can be enabled to ensure
that every piece has exactly two neighbours.
"""
function optim(
    pieces::Union{Board,Vector{Polymino}},
    omega::CellSet;
    opt::Dict = Dict(),
    highs::Bool = false,
    connect2::Bool = false,
    mask::Dict = Dict(),
)::Board

    model = highs ? Model(HiGHS.Optimizer) : Model(Gurobi.Optimizer)

    if !highs
        set_optimizer_attribute(model, "Threads", 1)
    end
    for (k, v) in opt
        @info "set $(k)=$(v) for optimizer"
        set_optimizer_attribute(model, k, v)
    end

    if pieces isa Board
        fixed_pieces = Set(p for p in pieces if !issubset(p.cells, omega))
        movable_pieces = Set(canonical.(setdiff(pieces, fixed_pieces)))
        xdom = union(omega, holes(pieces), cells(pieces))
    else
        fixed_pieces = Set{Polymino}()
        movable_pieces = Set(canonical.(pieces))
        xdom = union(omega, holes(omega))
    end

    poly = vcat(fixed_pieces..., movable_pieces...)

    idx = Dict(x => i for (i, x) in enumerate(xdom))


    polys = Polymino[]

    form = collect(Vector{Int}() for _ = 1:length(poly))
    unique = [p for p in pieces if count(==(canonical(p)), canonical.(pieces)) == 1]
    not_rotated = isempty(unique) ? Polymino(Pt[]) : argmax(p -> length(sym(p)), unique)
    @info "not rotated" not_rotated
    for (i, p) in enumerate(poly)
        if p in fixed_pieces
            push!(polys, p)
            push!(form[i], length(polys))
        else
            cur = Set(pp for xy in omega for pp in place_at(xy, p, omega, p != not_rotated))
            for pp in cur
                push!(polys, pp)
                push!(form[i], length(polys))
            end
        end
    end
    @info "nontrivial choices $(filter(x->x!=1, length.(form)))"

    np = length(polys)
    attached = Dict{Pt,Vector{Int}}(xy => Int[] for xy in xdom)
    for (i, pl) in enumerate(polys)
        for c in pl.cells
            push!(attached[c], i)
        end
    end

    ss = zeros(Int, np)
    xx = zeros(Int, length(xdom))
    @variable(model, s[1:np], Bin)
    @variable(model, x[1:length(xdom)], Bin)


    if pieces isa Board
        x0 = holes(pieces)
        for xy in xdom
            xx[idx[xy]] = xy in x0 ? 1 : 0
            set_start_value(x[idx[xy]], xy in x0 ? 1.0 : 0.0)
        end
        for (i, p) in enumerate(polys)
            if p in pieces
                set_start_value(s[i], 1.0)
                ss[i] = 1
            else
                set_start_value(s[i], 0.0)
            end
        end
    end

    for (i, p) in enumerate(poly)
        if p in fixed_pieces
            freqi = 1
        else
            freqi =
                count(==(canonical(p)), canonical.(pieces)) -
                count(==(canonical(p)), canonical.(fixed_pieces))
        end
        @constraint(model, sum(s[k] for k in form[i]) == freqi)
        if pieces isa Board && sum(ss[k] for k in form[i]) != freqi
            @info "violated" p sum(ss[k] for k in form[i]) freqi
        end
    end

    for xy in xdom
        @constraint(model, x[idx[xy]] + sum(s[k] for k in attached[xy]) <= 1)
        if pieces isa Board && xx[idx[xy]] + sum(ss[k] for k in attached[xy]; init = 0) > 1
            @info "violated" xy xx[idx[xy]] sum(ss[k] for k in attached[xy])
        end
        for nxy in n8(xy)
            if nxy in xdom
                @constraint(
                    model,
                    x[idx[xy]] <= x[idx[nxy]] + sum(s[k] for k in attached[nxy])
                )
            else
                @constraint(model, x[idx[xy]] == 0)
            end
        end
    end

    if connect2
        progress = 0
        for (i, p) in enumerate(polys)
            next_to_p = Set(vcat([attached[l] for l in hull(p) if haskey(attached, l)]...))
            # @constraint(model, s[i] => {2 == sum(s[j] for j in next_to_p)})
            @constraint(model, 10*s[i] + sum(s[j] for j in next_to_p)<=12)
        end
    end

    for l in xdom
        if haskey(mask, l)
            value = Int(mask[l])
            @constraint(model, x[idx[l]] == value)
        end
    end
 
    @objective(model, Max, sum(x[idx[xy]] for xy in xdom))

    optimize!(model)

    board = Board([polys[i] for i = 1:np if value(s[i]) > 0.5])
    nh = nholes(board)
    @info "after optim: nholes = $(nh)"
    board
end

"""
DFS to find an order in the board. If the pieces are already
on a ring, that will be detected and returned.
"""
function findchain(board::Board)::Vector{Polymino}
    queue = copy(board.pieces)
    chain = []
    while length(queue) > 0
        p = pop!(queue)
        if !(p in chain)
            push!(chain, p)
        end
        for q in board.pieces
            if any(c->c in p.cells, [c for x in q.cells for c in n4(x)]) &&
               q != p &&
               ! (q in chain)
                push!(chain, q)
                push!(queue, q)
            end
        end
    end
    chain
end

"""
Random processes that relaxes parts
of the current board and optimises only there.    
"""
function stochopt(
    board;
    chains = [2, 2, 2, 2, 2],
    stop = nothing,
    relax_by = 2,
    niter::Int = 500,
    opt::Dict = Dict(),
)::Board
    println(nholes(board))
    omega = relax(cells(board))
    xdom = union(omega, holes(board))
    best = nholes(board)
    np = length(board.pieces)
    for iter = 1:niter
        chain = findchain(board)
        x0 = holes(board)
        ch(sz) = begin
            i0 = rand(1:np)
            [chain[mod(i0+i, np)+1] for i = 1:sz]
        end
        mobile = vcat([ch(k) for k in chains]...)
        omega = relax(cells(Board(mobile)); relax_by = relax_by)
        board = optim(board, omega; opt = opt)
        @info "holes = $(nholes(board))"
        nh = nholes(board)
        if nh > best
            best = nh
            println(board)
        end
        if nh == stop
            @info "stop" nh iter
            break
        end
    end
    board
end
end
