#=
    This file includes all the funtions needed to compute the techno-economic
    indicators
 =#

#Power subscription prices interpatated/extrapolated to be on a continous scale.
# Prices comes from https://www.kelwatt.fr/guide/heures-creuses#:~:text=La%20p%C3%A9riode%20des%20heures%20creuses,et%20entre%2020h%20et%208h.
power_sub = [6, 9, 12, 15, 18, 24, 30, 36]
prices = [144.32, 183.63, 221.81, 258.98, 294.25, 360.61, 428.27, 494.92]
interp_linear_extrap_sub_prices = LinearInterpolation(power_sub, prices, extrapolation_bc=Line())

mutable struct COST{T <: Array{Float64}}
    capex::T
    opex::T
    total::T
end

# Compute costs
COST(mg::Microgrid, designer::AbstractDesigner) = COST(1:mg.parameters.ns, mg, designer)
# Compute costs for a given scenario s
function COST(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner)

    γ = repeat(1. ./ (1. + mg.parameters.τ) .^ range(0, length = mg.parameters.ny, step = mg.parameters.Δy), 1, length(s))

    capexx = γ .* capex(s, mg, designer)

    # Discounted opex
    opex = γ .* dropdims(complex_grid_cost(s, mg, designer),dims=1)

    return COST(capexx, opex, capexx .+ opex )
end

 mutable struct NPV{T <: Array{Float64}}
      capex::T
      opex::T
      salvage::T
      cf::T
      cumulative_npv::T
      total::T
 end

 # Compute costs
NPV(mg::Microgrid, designer::AbstractDesigner) = NPV(1:mg.parameters.ns, mg, designer)
# Compute costs for a given scenario s
function NPV(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner)

    # Discount factor
    γ = repeat(1. ./ (1. + mg.parameters.τ) .^ range(0, length = mg.parameters.ny, step = mg.parameters.Δy), 1, length(s))

    # Discounted capex
    capexx = γ .* capex(s, mg, designer)

    # Discounted opex
    #opex = γ .* dropdims(baseline_cost(s, mg) .-  grid_cost(s, mg, designer), dims=1)
    opex = γ .* dropdims(baseline_cost(s, mg) .-  complex_grid_cost(s, mg, designer), dims=1)

    # Discounted salvage value
    salvage = γ .*  salvage_value(s, mg)

    # Discounted cash flow
    cf = - capexx .+ opex .+ salvage

    # Discounted NPV each year
    cumulative = cumsum(cf, dims = 1)

    # Discounted NPV
    total = sum(cf, dims=1)

    return NPV(capexx, opex, salvage, cf, cumulative, total)
end

# Baseline cost
baseline_cost(mg::Microgrid) = baseline_cost(1:mg.parameters.ns, mg)
baseline_cost(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid) = baseline_cost(1:mg.parameters.ny, s, mg)
function baseline_cost(y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid)
    # TODO : compute baseline cost whatever the microgrid...



    total = 0.
    if !isempty(mg.demands)
        for (k,a) in enumerate(mg.demands)
            if a.carrier isa Electricity
                #
                total = total .+ sum(a.carrier.power[:,y,s] .* mg.grids[1].cost_in[:,y,s] * mg.parameters.Δh, dims = 1)

                #Hupothese l'abonnement était designé parfaitement pour l'utilisation
                # Entrainant 0 coût de depassement
                P_max = maximum(a.carrier.power[:,y,s], dims = 1)
                subscribtion = zeros(1,length(y), length(s))
                for i in y
                    for j in s
                        subscribtion[:,i,j] = interp_linear_extrap_sub_prices(P_max[:,i,j])
                    end
                end
                total = total .+ subscribtion  #abonnement


            elseif a.carrier isa Heat
                total = total .+ sum(a.carrier.power[:,y,s] / mg.converters[isin(mg.converters, Heater)[2]].η_E_H .* mg.grids[1].cost_in[:,y,s] * mg.parameters.Δh, dims = 1)
            end


        end
    end


    return total
end
# Grid cost
grid_cost(mg::Microgrid, designer::AbstractDesigner) = grid_cost(1:mg.parameters.ns, mg , designer)
grid_cost(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner) = grid_cost(1:mg.parameters.ny, s, mg, designer)
function grid_cost(years::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner)

    return sum(sum(max.(0., a.carrier.power[:,years,s]) .* a.cost_in[:,years,s] .- min.(0., a.carrier.power[:,years,s]) .* a.cost_out[:,years,s], dims = 1) * mg.parameters.Δh for a in mg.grids) # Energy buying cost

end


complex_grid_cost(mg::Microgrid, designer::AbstractDesigner) = complex_grid_cost(1:mg.parameters.ns, mg , designer)
complex_grid_cost(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner) = complex_grid_cost(1:mg.parameters.ny, s, mg, designer)
function complex_grid_cost(years::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner)

    hour_factor = ones(size(mg.grids[1].cost_in))
    creuse = (1:8760)[((1:8760).%24 .- 6 .< 0) .| ((1:8760).%24 .- 22 .>= 0)]
    hour_factor[creuse,:,:] .* 0.75   #heure pleine . #Informations comes from https://www.kelwatt.fr/guide/heures-creuses#:~:text=La%20p%C3%A9riode%20des%20heures%20creuses,et%20entre%2020h%20et%208h.

    #Energy buying cost
    net_energy_cost = sum(sum(max.(0., a.carrier.power[:,years,s]) .* a.cost_in[:,years,s] .* hour_factor[:,years,s] .- min.(0., a.carrier.power[:,years,s]) .* a.cost_out[:,years,s], dims = 1) * mg.parameters.Δh for a in mg.grids)

    #overcome cost by year
    overcome_cost = zeros(length(years),length(s))
    for s in 1:length(s)
        for y in 1:length(years)
            overcome_cost[y,s] = sum(sum(count(nb_overcome->(nb_overcome > 0), a.carrier.power[:,y,s] .- designer.decisions.subscribed_power[k][y,s]) * a.cost_exceed[y,s] for (k, a) in enumerate(mg.grids))) # count the hourly consumption exceeding the subscribed power
        end
    end

    net_energy_cost[1,:,:] = net_energy_cost[1,:,:] .+ overcome_cost
    return net_energy_cost
end

# CAPEX
capex(mg::Microgrid, designer::AbstractDesigner) = capex(1:mg.parameters.ns, mg, designer)
# CAPEX for a given scenario s
function capex(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner)



    capex = 0.
    # Generations
    for (k, a) in enumerate(mg.generations)
        capex = capex .+ designer.decisions.generations[k][:,s] .* a.cost[:,s]
    end
    # Storages
    for (k, a) in enumerate(mg.storages)
        capex = capex .+ designer.decisions.storages[k][:,s] .* a.cost[:,s]
    end
    # Converters
    for (k, a) in enumerate(mg.converters)
        capex = capex .+ designer.decisions.converters[k][:,s] .* a.cost[:,s]
    end


    # Subscribtion grid
    for (k, a) in enumerate(mg.grids)
        subscribtion = zeros(size(designer.decisions.subscribed_power[k][:,s]))
        for i in 1:size(designer.decisions.subscribed_power[k][:,s], 1)
            for j in s
                subscribtion[i,j] = interp_linear_extrap_sub_prices(designer.decisions.subscribed_power[k][i,j])
            end
        end
        capex = capex .+ subscribtion  #abonnement
    end
    return capex
end

# Salvage value
salvage_value(mg::Microgrid) = salvage_value(1:mg.parameters.ns, mg)
# Salvage value for a given scenario s
function salvage_value(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid)
    # Linear depreciation of components
    nh, ny = mg.parameters.nh, mg.parameters.ny
    salvage = zeros(mg.parameters.ny, length(s))
    salvage[ny,:] .= 0.
    # Generations
    for a in mg.generations
        salvage[ny,:] = salvage[ny,:] .+ (a.lifetime .- ny) ./ a.lifetime .* a.cost[ny, s] .* a.powerMax[ny,s]
    end
    # Storages
    for a in mg.storages
        if a isa AbstractLiion
            salvage[ny,:] = salvage[ny,:] .+ ((a.soh[1,end,s] .- a.SoH_threshold) ./ (1 .-a.SoH_threshold)) .* a.cost[ny, s] .* a.Erated[ny,s] # 100% value at 100% SOH, 0% at EOL
            #salvage[ny,:] = salvage[ny,:] .+ a.soh[1,end,s] .* a.cost[ny, s] .* a.Erated[ny,s]
            #$a.soh[end,end,s]$ remplace ici $(a.lifetime .- ny) ./ a.lifetime$ comme indicateur de la fraction de vie restante
        end
    end
    # Converters
    for a in mg.converters
        salvage[ny,:] = salvage[ny,:] .+ (a.lifetime .- ny) ./ a.lifetime .* a.cost[ny, s]
    end
    return salvage
end

mutable struct EAC{T <: Array{Float64}}
     capex::T
     opex::T
     total::T
end

function EAC(y::Int64, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner)
    # Annualised capex
    capex = Genesys.annualised_capex(1:y, s, mg, designer)
    # opex
    opex = grid_cost(y, s, mg, designer)

    return EAC(capex, opex, capex .+ opex)
end
# Annualised CAPEX
function annualised_capex(y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner)
    # Preallocation
    capex = 0.
    # Generations
    for (k, a) in enumerate(mg.generations)
        Γ = (mg.parameters.τ * (mg.parameters.τ + 1.) ^ a.lifetime) / ((mg.parameters.τ + 1.) ^ a.lifetime - 1.)
        capex = capex .+ Γ .* designer.decisions.generations[k][y,s] .* a.cost[y,s]
    end
    # Storages
    for (k, a) in enumerate(mg.storages)
        Γ = (mg.parameters.τ * (mg.parameters.τ + 1.) ^ a.lifetime) / ((mg.parameters.τ + 1.) ^ a.lifetime - 1.)
        capex = capex .+ Γ .* designer.decisions.storages[k][y,s] .* a.cost[y,s]
    end
    # Converters
    for (k, a) in enumerate(mg.converters)
        Γ = (mg.parameters.τ * (mg.parameters.τ + 1.) ^ a.lifetime) / ((mg.parameters.τ + 1.) ^ a.lifetime - 1.)
        capex = capex .+ Γ .* designer.decisions.converters[k][y,s] .* a.cost[y,s]
    end

    # Subscribtion grid
    for (k, a) in enumerate(mg.grids)
        tmp = zeros(size(designer.decisions.subscribed_power[k][:,s]))
        for i in 1:size(designer.decisions.subscribed_power[k][:,s],1)
            for j in s
                tmp[i,j] = interp_linear_extrap_sub_prices(designer.decisions.subscribed_power[k][i,j])
            end
        end
        capex = capex .+ (tmp/length(y))    #abonnement
    end
    return capex
end

# Share of renewables
renewable_share(mg::Microgrid) = renewable_share(1:mg.parameters.ns, mg)
# Share of renewables for a given scenario s
renewable_share(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid) = renewable_share(1:mg.parameters.ny, s, mg)
# Share of renewables for a given year y of a givn scenario s
function renewable_share(y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid)
    # TODO to be changed if there is no grid...
    total = 0.
    for (k,a) in enumerate(mg.demands)
        if a.carrier isa Electricity
            total = total .+ sum(a.carrier.power[:,y,s], dims = 1)
        elseif a.carrier isa Heat
            total = total .+ sum(a.carrier.power[:,y,s], dims = 1) ./ mg.converters[isin(mg.converters, Heater)[2]].η_E_H
        end
    end
    for (k,a) in enumerate(mg.grids)
        if a.carrier isa Electricity
            return share = dropdims(1. .- sum(max.(0., a.carrier.power[:,y,s]), dims = 1) ./ total, dims=1)
        else
            println("Renewable share not yet defined!")
            return nothing
        end
    end
end

# LPSP
mutable struct LPSP{T}
    elec::Union{Nothing, T}
    heat::Union{Nothing, T}
    EnergyCarrier::Union{Nothing, T}
end

LPSP(mg::Microgrid) = LPSP(1:mg.parameters.ns, mg)
# LPSP for a given scenario s
LPSP(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid) = LPSP(1:mg.parameters.ny, s, mg)
# LPSP for a given scenario s and year y
function LPSP(y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid)
    # Initialization
    elec, heat, EnergyCarrier = nothing, nothing, nothing
    # Computation
    for a in mg.demands
        if a.carrier isa Electricity
            elec = sum(max.(0., power_balance(1:mg.parameters.nh, y, s, mg, Electricity)), dims=1) ./ sum(a.carrier.power[:, y, s], dims = 1)
            for aa in mg.grids
                if aa.carrier isa Electricity
                    elec = elec .- sum(max.(0., aa.carrier.power[:, y, s]), dims=1) ./ sum(a.carrier.power[:, y, s], dims = 1)
                end
            end
            elec = dropdims(elec,dims=1)
        elseif a.carrier isa Heat
            heat = sum(max.(0., power_balance(1:mg.parameters.nh, y, s, mg, Heat)), dims=1) ./ sum(a.carrier.power[:, y, s], dims = 1)
            for aa in mg.grids
                if aa.carrier isa Heat
                    heat = heat .- sum(max.(0., aa.carrier.power[:, y, s]), dims=1) ./ sum(a.carrier.power[:, y, s], dims = 1)
                end
            end
            heat = dropdims(heat,dims=1)
        elseif a.carrier isa Hydrogen
            EnergyCarrier = sum(max.(0., power_balance(1:mg.parameters.nh, y, s, mg, Hydrogen)), dims=1) ./ sum(a.carrier.power[:, y, s], dims = 1)
            for aa in mg.grids
                if aa.carrier isa Electricity
                    hydrogen = hydrogen .- sum(max.(0., aa.carrier.power[:, y, s]), dims=1) ./ sum(a.carrier.power[:, y, s], dims = 1)
                end
            end
            hydrogen = dropdims(hydrogen,dims=1)
        end
    end
    return LPSP(elec, heat, EnergyCarrier)
end

mutable struct Metrics{T}
    baseline::T
    npv::NPV{T}
    eac::EAC{T}
    renewable_share::T
    lpsp::LPSP{T}
    cost::COST{T}
end

# Compute indicators
Metrics(mg::Microgrid, designer::AbstractDesigner) = Metrics(1:mg.parameters.ns, mg, designer)
# Compute indicators for a given scenario s
function Metrics(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner)
    # Baseline cost
    baseline = dropdims(baseline_cost(mg), dims = 1)
    # NPV
    npv = NPV(s, mg, designer)
    # EAC
    eac = EAC(mg.parameters.ny, s, mg, designer)
    # Share of renewables
    share = renewable_share(s, mg)
    # LPSP
    lpsp = LPSP(s, mg)

    cost = COST(s, mg, designer)

    return Metrics(baseline, npv, eac, share, lpsp, cost)
end
