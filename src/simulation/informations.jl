#=
    This file includes all the funtions needed to update the operation
    and investment informations
 =#

function update_operation_informations!(h::Int64, y::Int64, s::Int64, mg::Microgrid, ω::Scenarios)
    # Demands
    for (k, a) in enumerate(mg.demands)
        a.timestamp[h,y,s] = ω.demands[k].t[h,y,s]
        a.carrier.power[h,y,s] = ω.demands[k].power[h,y,s]
    end
    # Generations
    for (k, a) in enumerate(mg.generations)
        a.timestamp[h,y,s] = ω.generations[k].t[h,y,s]
        a.carrier.power[h,y,s] = a.powerMax[y,s] * ω.generations[k].power[h,y,s]
    end

end

function update_investment_informations!(y::Int64, s::Int64, mg::Microgrid, ω::Scenarios)
    # Generations
    for (k, a) in enumerate(mg.generations)
        a.cost[y,s] = ω.generations[k].cost[y,s]
    end
    # Storages
    for (k, a) in enumerate(mg.storages)
        a.cost[y,s] = ω.storages[k].cost[y,s]
    end
    # Converters
    for (k, a) in enumerate(mg.converters)
        a.cost[y,s] = ω.converters[k].cost[y,s]
    end

end


function update_grid_cost_informations!(y::Int64, s::Int64, mg::Microgrid, ω::Scenarios)
    # Grids - We assume the price of electricity is known over the year
    for (k, a) in enumerate(mg.grids)
            a.cost_in[1:end,y,s] = ω.grids[k].cost_in[1:end,y,s]
            a.cost_out[1:end,y,s] = ω.grids[k].cost_out[1:end,y,s]
            a.cost_exceed[y,s] = ω.grids[k].cost_exceed[y,s]
    end

end
