#=
    Loads modelling
 =#

mutable struct Demand <: AbstractDemand
     # Variables
     carrier::EnergyCarrier
     timestamp::Array{DateTime,3}

     # Inner constructor
     Demand(; carrier = Electricity()) = new(carrier)
end

### Preallocation
function preallocate!(ld::Demand, nh::Int64, ny::Int64, ns::Int64)
    ld.carrier.power = convert(SharedArray,zeros(nh, ny, ns))
    ld.timestamp = Array{DateTime}(undef,(nh, ny, ns))
    return ld
end
