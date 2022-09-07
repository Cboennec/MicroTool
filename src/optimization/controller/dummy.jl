#=
    Dummy controller
=#

mutable struct Dummy <: AbstractController
    decisions::NamedTuple
    Dummy() = new()
end

### Offline
function initialize_controller!(mg::Microgrid, controller::Dummy, ω::AbstractScenarios)
    # Preallocation
    preallocate!(mg, controller)
    return controller
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::Dummy)
    return controller
end
