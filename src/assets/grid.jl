#=
    Elec grid modelling
 =#

mutable struct Grid <: AbstractGrid
    # Parameters
    powerMax::Float64
    # Variables
    carrier::EnergyCarrier
    cost_in::AbstractArray{Float64,3}
    cost_out::AbstractArray{Float64,3}
    cost_exceed::AbstractArray{Float64,2}

    # Inner constructor
    Grid(; powerMax = 36., carrier = Electricity()) = new(powerMax, carrier)
end

### Preallocation
function preallocate!(grid::Grid, nh::Int64, ny::Int64, ns::Int64)
     grid.carrier.power = convert(SharedArray,zeros(nh, ny, ns))
     grid.cost_in = convert(SharedArray,zeros(nh, ny, ns))
     grid.cost_out = convert(SharedArray,zeros(nh, ny, ns))
     grid.cost_exceed = convert(SharedArray, ones(ny, ns) .* 11. ) # price from https://electricitedesavoie.fr/2017/05/09/3-points-comprendre-eviter-depassement-de-puissance-souscrite/#:~:text=Cette%20puissance%20est%20exprim%C3%A9e%20en,%2C11%E2%82%AC%2Fheure%20d%C3%A9pass%C3%A9e.
     return grid
end
