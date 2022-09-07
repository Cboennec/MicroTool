#=
    This file includes all the functions needed to check the power balance
    constraints
=#


#TODO Here is a sub optimal implementation with the types comparison it's too slow there must be a way to write a bit more code that go really faster.

function compute_power_balances!(h::Int64, y::Int64, s::Int64, mg::Microgrid)
    # Hydrogen
    checking!(h, y, s, mg, Hydrogen)
    # Heat
    checking!(h, y, s, mg, Heat)
    # Electricity
    checking!(h, y, s, mg, Electricity)
end

function power_balance(h::Union{Int64, UnitRange{Int64}}, y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, type::DataType)
    # Energy balance
    balance = 0.
    # Demands
    if !isempty(mg.demands)
        for a in mg.demands
            a.carrier isa type ? balance = balance .+ a.carrier.power[h,y,s] : nothing
        end
    end
    # Generations
    if !isempty(mg.generations)
        for a in mg.generations
            a.carrier isa type ? balance = balance .- a.carrier.power[h,y,s] : nothing
        end
    end
    # Storages
    if !isempty(mg.storages)
        for a in mg.storages
            a.carrier isa type ? balance = balance .- a.carrier.power[h,y,s] : nothing
        end
    end
    # Converters
    if !isempty(mg.converters)
        for a in mg.converters
            for c in a.carrier
                c isa type ? balance = balance .- c.power[h,y,s] : nothing
            end
        end
    end
    return balance
end

function power_balance(h::Int64, y::Int64, s::Int64, mg::Microgrid, type::DataType)
    # Energy balance
    balance = 0.
    # Demands
    if !isempty(mg.demands)
        for a in mg.demands
            a.carrier isa type ? balance = balance + a.carrier.power[h,y,s] : nothing
        end
    end
    # Generations
    if !isempty(mg.generations)
        for a in mg.generations
            a.carrier isa type ? balance = balance - a.carrier.power[h,y,s] : nothing
        end
    end
    # Storages
    if !isempty(mg.storages)
        for a in mg.storages
            a.carrier isa type ? balance = balance - a.carrier.power[h,y,s] : nothing
        end
    end
    # Converters
    if !isempty(mg.converters)
        for a in mg.converters
            for c in a.carrier
                c isa type ? balance = balance - c.power[h,y,s] : nothing
            end
        end
    end
    return balance
end

function checking!(h::Int64, y::Int64, s::Int64, mg::Microgrid, type::DataType)
    # Energy balance
    balance = power_balance(h, y, s, mg, type)
    # Grids
    if !isempty(mg.grids)
        for a in mg.grids
            if a.carrier isa type
                a.carrier.power[h,y,s] = min(a.powerMax, max(-a.powerMax, balance))
            end
        end
    else # If the balance equation is not fulfilled, systems are turned to zerp
        if balance > sum(a.carrier.power[h,y,s] for a in mg.demands if a.carrier isa type) + 1e-2 # 0.01 kW tolerance
            # Storage set to zero
            for a in mg.storages
                if a.carrier isa type
                    a.carrier.power[h,y,s]= 0.
                    a.soc[h+1,y,s] = max(0., a.soc[h,y,s] * (1. - a.η_self))
                end
            end
            # Converters set to zero
            for a in mg.converters
                for c in a.carrier
                    c.power[h,y,s] = 0.
                end
            end
        else
            println("Shedding energy demand!")
        end
    end
end
