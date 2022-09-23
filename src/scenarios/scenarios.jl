#=
    Scenario reduction functions
=#
abstract type AbstractScenarios end

mutable struct Scenarios{T, O, I} <: AbstractScenarios
    demands::Vector{NamedTuple{(:t, :power),Tuple{T,O}}}
    generations::Vector{NamedTuple{(:t, :power, :cost), Tuple{T, O, I}}}
    storages::Vector{NamedTuple{(:cost,), Tuple{I}}}
    converters::Vector{NamedTuple{(:cost,), Tuple{I}}}
    grids::Vector{NamedTuple{(:cost_in, :cost_out, :cost_exceed), Tuple{O, O, I}}}
end

# Constructor
function Scenarios(mg::Microgrid, d::Dict{})
    # Utils to simplify the writting
    h, y, s = 1:mg.parameters.nh, 1:mg.parameters.ny, 1:mg.parameters.ns
    T, O, I = Array{DateTime,3}, Array{Float64, 3}, Array{Float64, 2}

    # Initialize
    demands = Vector{NamedTuple{(:t, :power),Tuple{T,O}}}(undef, length(mg.demands))
    generations = Vector{NamedTuple{(:t, :power, :cost), Tuple{T, O, I}}}(undef, length(mg.generations))
    storages = Vector{NamedTuple{(:cost,), Tuple{I}}}(undef, length(mg.storages))
    converters = Vector{NamedTuple{(:cost,), Tuple{I}}}(undef, length(mg.converters))
    grids = Vector{NamedTuple{(:cost_in, :cost_out, :cost_exceed), Tuple{O, O, I}}}(undef, length(mg.grids))
    # Demands
    for (k, a) in enumerate(mg.demands)
        if a.carrier isa Electricity
            demands[k] = (t = d["ld_E"].t[h, y, s], power = d["ld_E"].power[h, y, s])
        elseif a.carrier isa Heat
            demands[k] = (t = d["ld_H"].t[h, y, s], power = d["ld_H"].power[h, y, s])
        end
    end
    # Generation
    for (k, a) in enumerate(mg.generations)
        if a isa Solar
            generations[k] = (t = d["pv"].t[h, y, s], power = d["pv"].power[h, y, s], cost = d["pv"].cost[y, s])
        end
    end
    # Storages
    for (k, a) in enumerate(mg.storages)

        if typeof(a) <: AbstractLiion
            storages[k] = (cost = d["liion"].cost[y, s],)
        elseif a isa ThermalStorage
            storages[k] = (cost = d["tes"].cost[y, s],)
        elseif a isa H2Tank
            storages[k] = (cost = d["h2tank"].cost[y, s],)
        end
    end
    # Converters
    for (k, a) in enumerate(mg.converters)
        if a isa Electrolyzer
            converters[k] = (cost = d["elyz"].cost[y, s],)
        elseif a isa FuelCell
            converters[k] = (cost = d["fc"].cost[y, s],)
        elseif a isa Heater
            converters[k] = (cost = d["heater"].cost[y, s],)
        end
    end
    # Grids
    for (k, a) in enumerate(mg.grids)
        if a.carrier isa Electricity
            grids[k] = (cost_in = d["grid"].cost_in[h, y, s], cost_out = d["grid"].cost_out[h, y, s], cost_exceed = zeros(length(y),length(s)) .+ 10) #TODO this price should come from the scenarios
        end
    end

    return Scenarios(demands, generations, storages, converters, grids)
end


# repetitive year for longer scenarios
function Scenarios(mg::Microgrid, d::Dict{}; same_year = false, seed = []) # repeat make every year the same, seed decide with year to use.
    # Utils to simplify the writting


    h, y, s = 1:mg.parameters.nh, 2:2, 1:mg.parameters.ns
    T, O, I = Array{DateTime,3}, Array{Float64, 3}, Array{Float64, 2}

    rep_time = convert(Int64,  mg.parameters.ny)
    # Initialize
    demands = Vector{NamedTuple{(:t, :power),Tuple{T,O}}}(undef, length(mg.demands))
    generations = Vector{NamedTuple{(:t, :power, :cost), Tuple{T, O, I}}}(undef, length(mg.generations))
    storages = Vector{NamedTuple{(:cost,), Tuple{I}}}(undef, length(mg.storages))
    converters = Vector{NamedTuple{(:cost,), Tuple{I}}}(undef, length(mg.converters))
    grids = Vector{NamedTuple{(:cost_in, :cost_out, :cost_exceed), Tuple{O, O, I}}}(undef, length(mg.grids))
    # Demands
    for (k, a) in enumerate(mg.demands)
        if a.carrier isa Electricity
            demands[k] = (t = repeat(d["ld_E"].t[h, y, s], outer = (1,rep_time,1)), power = compose(d["ld_E"].power, rep_time, mg, 3; rep = same_year, s_num = seed))
        elseif a.carrier isa Heat
            demands[k] = (t = repeat(d["ld_H"].t[h, y, s], outer = (1,rep_time,1)), power = compose(d["ld_H"].power, rep_time, mg, 3; rep = same_year, s_num = seed))
        end
    end
    # Generation
    for (k, a) in enumerate(mg.generations)
        if a isa Solar
            #generations[k] = (t = d["pv"].t[h, y, s], power = d["pv"].power[h, y, s], cost = d["pv"].cost[y, s])
            generations[k] = (t = repeat(d["pv"].t[h, y, s], outer = (1,rep_time,1)), power = compose( d["pv"].power, rep_time, mg, 3; rep = same_year, s_num = seed), cost = compose(d["pv"].cost, rep_time, mg, 2;  rep = same_year, s_num = seed))
        end
    end
    # Storages
    for (k, a) in enumerate(mg.storages)

        if typeof(a) <: AbstractLiion
        #    storages[k] = (cost = d["liion"].cost[y, s],)
             storages[k] = (cost = compose(d["liion"].cost, rep_time, mg, 2; rep = same_year, s_num = seed),)
        elseif a isa ThermalStorage
        #    storages[k] = (cost = d["tes"].cost[y, s],)
            storages[k] = (cost = compose(d["tes"].cost, rep_time, mg, 2; rep = same_year, s_num = seed),)
        elseif a isa H2Tank
        #    storages[k] = (cost = d["h2tank"].cost[y, s],)
            storages[k] = (cost = compose(d["h2tank"].cost, rep_time, mg, 2; rep = same_year, s_num = seed),)
        end
    end
    # Converters
    for (k, a) in enumerate(mg.converters)
        if a isa Electrolyzer
            #converters[k] = (cost = d["elyz"].cost[y, s],)
            converters[k] = (cost = compose(d["elyz"].cost, rep_time, mg, 2; rep = same_year, s_num = seed))
        elseif a isa FuelCell
            #converters[k] = (cost = d["fc"].cost[y, s],)
            converters[k] = (cost = compose(d["fc"].cost, rep_time, mg, 2; rep = same_year, s_num = seed))

        elseif a isa Heater
            #converters[k] = (cost = d["heater"].cost[y, s],)
            converters[k] = (cost = compose(d["heater"].cost, rep_time, mg, 2; rep = same_year, s_num = seed))
        end
    end
    # Grids
    for (k, a) in enumerate(mg.grids)
        if a.carrier isa Electricity
            #grids[k] = (cost_in = d["grid"].cost_in[h, y, s], cost_out = d["grid"].cost_out[h, y, s])
            grids[k] = (cost_in = compose(d["grid"].cost_in, rep_time, mg, 3; rep = same_year, s_num = seed), cost_out = compose(d["grid"].cost_out, rep_time, mg, 3; rep = same_year, s_num = seed), cost_exceed = zeros( mg.parameters.ny,  mg.parameters.ns) .+ 10)#TODO this price should come from the scenarios
        end
    end

    return Scenarios(demands, generations, storages, converters, grids)
end


#Compose the data for a longer scenario based on existing scenario of 1 year
#The dim parameter defines the number of dimension on which this data is represented
#optionnal
#The rep parameter defines wether or not every year have to be the same
#s_num defines the id of the scenarios we are going to use, very useful for reproductivity of the results
# If no seed is provided, the scenario IDs will  be randomly selected.
function compose(array, rep_time::Int64, mg::Microgrid, dim::Int64; rep = false, s_num = [] )
    nh = mg.parameters.nh
    ny = mg.parameters.ny
    ns = mg.parameters.ns

    hours = 1:nh
    years = 1:ny
    scenarios = 1:ns

    r = 0

    if dim == 3
        result = repeat(array[hours, 2:2, scenarios], outer = (1,rep_time,1)) # instantiate an array of the right size.

        for s in scenarios
            if !isempty(s_num)
                r =  s_num[s]
            else
                r = convert(Int64, floor(rand() * 1000)+1)
            end

            for y in years
                if !rep #On refait un tirage
                    r = convert(Int64, floor(rand() * 1000)+1)
                end
                result[:,y,s] = array[:,2,r] # Get the  whole second year of a random scenario and plug in the selected year
            end
        end
    elseif dim == 2
        result = repeat(array[2:2, scenarios], outer = (rep_time,1))

        for s in scenarios
            if !isempty(s_num)
                r =  s_num[s]
            else
                r = convert(Int64, floor(rand() * 1000)+1)
            end

            for y in years
                if !rep #On refait un tirage
                    r = convert(Int64, floor(rand() * 1000)+1)
                end
                result[y,s] = array[2,r] #Get the  whole second year of a random scenario
            end
        end
    end

    return result

end
