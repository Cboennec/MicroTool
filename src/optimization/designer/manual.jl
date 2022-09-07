#=
    Manual designer
=#

mutable struct Manual <: AbstractDesigner
    generations::Vector{Float64}
    storages::Vector{Float64}
    converters::Vector{Float64}
    subscribed_power::Vector{Float64}
    decisions::NamedTuple

    Manual(; generations = [0.], storages = [0.], converters = [0.], subscribed_power = [0.]) = new(generations, storages, converters, subscribed_power)
end

### Offline
function initialize_designer!(mg::Microgrid, designer::Manual, ω::AbstractScenarios)
    # Preallocation
    preallocate!(mg, designer)

    # Fix initial values

    for k in 1:length(mg.generations)
        designer.decisions.generations[k][1,:] .= designer.generations[k]
    end
    for k in 1:length(mg.storages)
        designer.decisions.storages[k][1,:] .=  designer.storages[k]
    end
    for k in 1:length(mg.converters)
        designer.decisions.converters[k][1,:] .= designer.converters[k]
    end
    for k in 1:length(mg.grids)
        designer.decisions.subscribed_power[k][:,:] .= designer.subscribed_power[k]
    end

    return designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, mg::Microgrid, designer::Manual)

    if mg.storages[1].soh[end,y,s] <= mg.storages[1].SoH_threshold
        designer.decisions.storages[1][y,s] = designer.storages[1]
    end
    return nothing
end
