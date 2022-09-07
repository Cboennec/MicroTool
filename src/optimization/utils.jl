# Risk measures
abstract type AbstractRiskMeasure end
# Expectation
struct Expectation <: AbstractRiskMeasure end
# Worst case
struct WorstCase <: AbstractRiskMeasure end
# Conditional value at risk
struct CVaR <: AbstractRiskMeasure
    β::Float64

    function CVaR(β::Float64)
        @assert 0. <= β <= 1. "β must be in [0,1]"
        return new(β)
    end
end

conditional_value_at_risk(support::Array{Float64,1}, probabilities::Array{Float64,1}, risk::WorstCase) = conditional_value_at_risk(support, probabilities, 0.)
conditional_value_at_risk(support::Array{Float64,1}, probabilities::Array{Float64,1}, risk::Expectation) = conditional_value_at_risk(support, probabilities, 1.)
conditional_value_at_risk(support::Array{Float64,1}, probabilities::Array{Float64,1}, risk::CVaR) = conditional_value_at_risk(support, probabilities, 1. - risk.β)

# From https://github.com/jaantollander/ConditionalValueAtRisk
function conditional_value_at_risk(support::Array{Float64,1}, probabilities::Array{Float64,1}, α::Float64)
    # Value at risk
    var = value_at_risk(support, probabilities, α)
    # Tail values
    if α == 0.
        return var
    else
        tail = support .< var
        return (sum(probabilities[tail] .* support[tail]) - (sum(probabilities[tail]) - α) * var) / α
    end
end

function value_at_risk(support::Array{Float64,1}, probabilities::Array{Float64,1}, α::Float64)

    i = findfirst(cumsum(probabilities[sortperm(support)]) .>= α)
    if i == nothing
        return sort(support)[end]
    else
        return sort(support)[i]
    end
end

# MILP functions
# Decisions
#Investment on the energy generation units (define the size)
function add_investment_decisions!(m::Model, generations::Vector{AbstractGeneration})
    if !isempty(generations)
        na = length(generations)
        @variable(m, r_g[1:na])
    end
end
#Investment on the energy storage units (define the size)
function add_investment_decisions!(m::Model, storages::Vector{AbstractStorage})
    if !isempty(storages)
        na = length(storages)
        @variable(m, r_sto[1:na])
    end
end
#Investment on the energy converter units
function add_investment_decisions!(m::Model, converters::Vector{AbstractConverter})
    if !isempty(converters)
        na = length(converters)
        @variable(m, r_c[1:na])
    end
end

#Multi year dynamic version
function add_investment_decisions!(m::Model, generations::Vector{AbstractGeneration}, ny::Integer)
    if !isempty(generations)
        na = length(generations)
        @variable(m, r_g[1:na])
    end
end
function add_investment_decisions!(m::Model, storages::Vector{AbstractStorage}, ny::Integer)
    if !isempty(storages)
        na = length(storages)
        @variable(m, r_sto[1:ny, 1:na])
    end
end
function add_investment_decisions!(m::Model, converters::Vector{AbstractConverter}, ny::Integer)
    if !isempty(converters)
        na = length(converters)
        @variable(m, r_c[1:ny, 1:na])
    end
end


function fix_investment_decisions!(m::Model, generations::Vector{Float64}, storages::Vector{Float64}, converters::Vector{Float64}, mg::Microgrid)
    # Generation
    if !isempty(mg.generations)
        fix.(m[:r_g], generations)
    end
    # Storages
    if !isempty(mg.storages)
        fix.(m[:r_sto], storages)
    end
    # Converters
    if !isempty(mg.converters)
        fix.(m[:r_c], converters)
    end
end

function fix_investment_decisions!(m::Model, generations::Vector{Float64}, storages::Vector{Float64}, converters::Vector{Float64})
    # Generation
    if !isempty(generations)
        fix.(m[:r_g], generations)
    end
    # Storages
    if !isempty(storages)
        fix.(m[:r_sto], storages)
    end
    # Converters
    if !isempty(converters)
        fix.(m[:r_c], converters)
    end
end


function add_operation_decisions!(m::Model, demands::Vector{AbstractDemand}, nh::Int64, ns::Int64)
    if !isempty(demands)
        na = length(demands)
        @variables(m, begin
        p_d[1:nh, 1:ns, 1:na]
        end)
    end
end
function add_operation_decisions!(m::Model, generations::Vector{AbstractGeneration}, nh::Int64, ns::Int64)
    if !isempty(generations)
        na = length(generations)
        @variables(m, begin
        p_g[1:nh, 1:ns, 1:na]
        end)
    end
end
function add_operation_decisions!(m::Model, storages::Vector{AbstractStorage}, nh::Int64, ns::Int64)
    if !isempty(storages)
        na = length(storages)
        @variables(m, begin
        p_ch[1:nh, 1:ns, 1:na]   >= 0.
        p_dch[1:nh, 1:ns, 1:na]  >= 0.
        soc[1:nh+1, 1:ns, 1:na]
        end)
    end
end
function add_operation_decisions!(m::Model, converters::Vector{AbstractConverter}, nh::Int64, ns::Int64)
    if !isempty(converters)
        na = length(converters)
        @variable(m, p_c[1:nh, 1:ns, 1:na] >= 0.)
    end
end
function add_operation_decisions!(m::Model, grids::Vector{AbstractGrid}, nh::Int64, ns::Int64)
    if !isempty(grids)
        na = length(grids)
        @variables(m, begin
        p_in[1:nh, 1:ns, 1:na]   >= 0.
        p_out[1:nh, 1:ns, 1:na]  >= 0.
        end)
    end
end

# Multi year dynamic version
function add_operation_decisions!(m::Model, demands::Vector{AbstractDemand}, nh::Int64, ns::Int64, ny::Int64)
    if !isempty(demands)
        na = length(demands)
        @variables(m, begin
        p_d[1:nh, 1:ns, 1:na, 1:ny]
        end)
    end
end
function add_operation_decisions!(m::Model, generations::Vector{AbstractGeneration}, nh::Int64, ns::Int64, ny::Int64)
    if !isempty(generations)
        na = length(generations)
        @variables(m, begin
        p_g[1:nh, 1:ny, 1:ns, 1:na]
        end)
    end
end
function add_operation_decisions!(m::Model, storages::Vector{AbstractStorage}, nh::Int64, ns::Int64, ny::Int64)
    if !isempty(storages)
        na = length(storages)
        @variables(m, begin
        p_ch[1:nh, 1:ny, 1:ns, 1:na]   >= 0.
        p_dch[1:nh, 1:ny, 1:ns, 1:na]  >= 0.
        soc[1:nh+1, 1:ny, 1:ns, 1:na]  >= 0.
        soh[1:nh+1, 1:ny, 1:ns, 1:na]  >= 0. #could exist over multiple storage here is just a battery
        end)
    end
end
function add_operation_decisions!(m::Model, converters::Vector{AbstractConverter}, nh::Int64, ns::Int64, ny::Int64)
    if !isempty(converters)
        na = length(converters)
        @variable(m, p_c[1:nh, 1:ny, 1:ns, 1:na] >= 0.)
    end
end
function add_operation_decisions!(m::Model, grids::Vector{AbstractGrid}, nh::Int64, ns::Int64, ny::Int64)
    if !isempty(grids)
        na = length(grids)
        @variables(m, begin
        p_in[1:nh, 1:ny, 1:ns, 1:na]   >= 0.
        p_out[1:nh, 1:ny, 1:ns, 1:na]  >= 0.
        end)
    end
end


# Multiyear case statz variables
function add_state_variables!(m::Model, ny::Int64)
    @variables(m, begin
    E_state[1:ny] >= 0.
    PV_state[1:ny] >= 0.
    end)
end


function add_invest_variables!(m::Model, ny::Int64)
    @variables(m, begin
    invest_B[1:ny], Bin
    #invest_PV[1:ny], Bin
    end)
end

# Investment bounds
function add_investment_constraints!(m::Model, generations::Vector{AbstractGeneration})
    if !isempty(generations)
        na = length(generations)
        @constraints(m, begin
        [a in 1:na], m[:r_g][a] >= generations[a].bounds.lb
        [a in 1:na], m[:r_g][a] <= generations[a].bounds.ub
        end)
    end
end
function add_investment_constraints!(m::Model, storages::Vector{AbstractStorage})
    if !isempty(storages)
        na = length(storages)
        @constraints(m, begin
        [a in 1:na], m[:r_sto][a] >= storages[a].bounds.lb
        [a in 1:na], m[:r_sto][a] <= storages[a].bounds.ub
        end)
    end
end
function add_investment_constraints!(m::Model, converters::Vector{AbstractConverter})
    if !isempty(converters)
        na = length(converters)
        @constraints(m, begin
        [a in 1:na], m[:r_c][a] >= converters[a].bounds.lb
        [a in 1:na], m[:r_c][a] <= converters[a].bounds.ub
        end)
    end
end

#Multi year dynamic version
function add_investment_constraints!(m::Model, generations::Vector{AbstractGeneration}, ny::Int64)
    if !isempty(generations)
        na = length(generations)
        @constraints(m, begin
        [a in 1:na], m[:r_g][a] >= generations[a].bounds.lb
        [a in 1:na], m[:r_g][a] <= generations[a].bounds.ub
        # [y in 1:ny, a in 1:na], m[:r_g][y,a] >= generations[a].bounds.lb
        # [y in 1:ny, a in 1:na], m[:r_g][y,a] <= generations[a].bounds.ub
        # [y in 1:ny, a in 1:na], m[:PV_state][y,a] <= generations[a].bounds.ub
        # [y in 1:ny, a in 1:na], m[:PV_state][y,a] >= generations[a].bounds.lb
        #
        # [y in 1:(ny-1), a in 1:na], m[:PV_state][y+1] - m[:r_g][y,a] >= -generations[a].bounds.ub * (1 - m[:invest_PV][y])
        # [y in 1:(ny-1), a in 1:na], m[:PV_state][y+1] - m[:r_g][y,a] <= generations[a].bounds.ub * (1 - m[:invest_PV][y])
        #
        # [y in 1:(ny-1), a in 1:na], m[:PV_state][y+1] - m[:PV_state][y] >= -generations[a].bounds.ub * m[:invest_PV][y]
        # [y in 1:(ny-1), a in 1:na], m[:PV_state][y+1] - m[:PV_state][y] <= generations[a].bounds.ub * m[:invest_PV][y]
        end)
    end
end
function add_investment_constraints!(m::Model, storages::Vector{AbstractStorage}, ny::Int64)
    if !isempty(storages)
        na = length(storages)
        @constraints(m, begin
        [y in 1:ny, a in 1:na], m[:r_sto][y,a] >= storages[a].bounds.lb
        [y in 1:ny, a in 1:na], m[:r_sto][y,a] <= storages[a].bounds.ub
        [y in 1:ny, a in 1:na], m[:E_state][y] <= storages[a].bounds.ub
        [y in 1:ny, a in 1:na], m[:E_state][y] >= storages[a].bounds.lb

        [y in 1:(ny-1), a in 1:na], m[:E_state][y+1] - m[:r_sto][y,a] >= -storages[a].bounds.ub * (1 - m[:invest_B][y])
        [y in 1:(ny-1), a in 1:na], m[:E_state][y+1] - m[:r_sto][y,a] <= storages[a].bounds.ub * (1 - m[:invest_B][y])

        [y in 1:(ny-1), a in 1:na], m[:E_state][y+1] - m[:E_state][y] >= -storages[a].bounds.ub * m[:invest_B][y]
        [y in 1:(ny-1), a in 1:na], m[:E_state][y+1] - m[:E_state][y] <= storages[a].bounds.ub * m[:invest_B][y]

        end)
    end
end
function add_investment_constraints!(m::Model, converters::Vector{AbstractConverter}, ny::Int64)
    if !isempty(converters)
        na = length(converters)
        @constraints(m, begin
        [y in 1:ny, a in 1:na], m[:r_c][y,a] >= converters[a].bounds.lb
        [y in 1:ny, a in 1:na], m[:r_c][y,a] <= converters[a].bounds.ub
        end)
    end
end



# Technical constraint
function add_technical_constraints!(m::Model, storages::Vector{AbstractStorage}, Δh::Int64, nh::Int64, ns::Int64)
    if !isempty(storages)
        na = length(storages)
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_dch][h,s,a] <= storages[a].α_p_dch * m[:r_sto][a]
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_ch][h,s,a]  <= storages[a].α_p_ch * m[:r_sto][a]
        # SoC bounds
        [h in 1:nh+1, s in 1:ns, a in 1:na], m[:soc][h,s,a] <= storages[a].α_soc_max * m[:r_sto][a]
        [h in 1:nh+1, s in 1:ns, a in 1:na], m[:soc][h,s,a] >= storages[a].α_soc_min * m[:r_sto][a]
        # State dynamics
        [h in 1:nh, s in 1:ns, a in 1:na], m[:soc][h+1,s,a] == m[:soc][h,s,a] * (1. - storages[a].η_self * Δh) - (m[:p_dch][h,s,a] / storages[a].η_dch - m[:p_ch][h,s,a] * storages[a].η_ch) * Δh
        # Initial and final states
        soc_ini[s in 1:ns, a in 1:na], m[:soc][1,s,a] == storages[a].soc_ini * m[:r_sto][a]
        end)
    end
end
function add_technical_constraints!(m::Model, converters::Vector{AbstractConverter}, nh::Int64, ns::Int64)
    if !isempty(converters)
        na = length(converters)
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_c][h,s,a]  <= m[:r_c][a]
        end)
    end
end
function add_technical_constraints!(m::Model, grids::Vector{AbstractGrid}, nh::Int64, ns::Int64)
    if !isempty(grids)
        na = length(grids)
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_in][h,s,a]  <= grids[a].powerMax
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_out][h,s,a] <= grids[a].powerMax
        end)
    end
end


# multi year dynamic version
function add_technical_constraints!(m::Model, storages::Vector{AbstractStorage}, Δh::Int64, nh::Int64, ny::Int64, ns::Int64)
    if !isempty(storages)
        na = length(storages)
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, y in 1:ny, s in 1:ns, a in 1:na], m[:p_dch][h,y,s,a] <= storages[a].α_p_dch * m[:r_sto][y,a]
        [h in 1:nh, y in 1:ny, s in 1:ns, a in 1:na], m[:p_ch][h,y,s,a]  <= storages[a].α_p_ch * m[:r_sto][y,a]
        # SoC bounds
        [h in 1:nh+1, y in 1:ny, s in 1:ns, a in 1:na], m[:soc][h,y,s,a] <= storages[a].α_soc_max * m[:r_sto][y,a]
        [h in 1:nh+1, y in 1:ny, s in 1:ns, a in 1:na], m[:soc][h,y,s,a] >= storages[a].α_soc_min * m[:r_sto][y,a]
        # State dynamics
        [h in 1:nh, y in 1:ny, s in 1:ns, a in 1:na], m[:soc][h+1,y,s,a] == m[:soc][h,y,s,a] * (1. - storages[a].η_self * Δh) - (m[:p_dch][h,y,s,a] / storages[a].η_dch - m[:p_ch][h,y,s,a] * storages[a].η_ch) * Δh
        # Initial and final states
        soc_ini[s in 1:ns, a in 1:na], m[:soc][1,1,s,a] == storages[a].soc_ini * m[:r_sto][a]

        #SoH evolution
        [h in 1:nh, y in 1:ny, s in 1:ns, a in 1:na], m[:soh][h+1,y,s] == m[:soh][h,y,s] -  (m[:p_dch][h,y,s,a] + m[:p_ch][h,y,s,a]) * Δh
        #SoH boundaries
        [h in 1:nh, y in 1:ny, s in 1:ns], m[:soh][h+1,y,s] >= (2 * storages[1].nCycle * m[:E_state][y]) * 0.9 #TODO remplacement de la batterie
        [h in 1:nh, y in 1:ny, s in 1:ns], m[:soh][h+1,y,s] <= 2 * storages[1].nCycle * m[:E_state][y]
        #If a new battery is bought it's assumed to be fully charge <ith a full SoH

        #Full batterie or last SoC
        [y in 1:(ny-1), s in 1:ns, a in 1:na], m[:soc][1,y+1,s] - (storages[a].α_soc_max * m[:r_sto][y,1]) <= storages[a].bounds.ub * storages[a].α_soc_max * (1-m[:invest_B][y])
        [y in 1:(ny-1), s in 1:ns, a in 1:na], m[:soc][1,y+1,s] - (storages[a].α_soc_max * m[:r_sto][y,1]) >= -storages[a].bounds.ub * storages[a].α_soc_max * (1-m[:invest_B][y])
        [y in 1:(ny-1), s in 1:ns, a in 1:na], m[:soc][1,y+1,s] - m[:soc][end,y,s] <= storages[a].bounds.ub * storages[a].α_soc_max * m[:invest_B][y]
        [y in 1:(ny-1), s in 1:ns, a in 1:na], m[:soc][1,y+1,s] - m[:soc][end,y,s] >= -storages[a].bounds.ub * storages[a].α_soc_max * m[:invest_B][y]
        #Full SOH or last SOH value
        [y in 1:(ny-1), s in 1:ns, a in 1:na], m[:soh][1,y+1,s] - 2 * storages[a].nCycle * m[:E_state][y] <= storages[a].bounds.ub * storages[a].α_soc_max * (1-m[:invest_B][y])
        [y in 1:(ny-1), s in 1:ns, a in 1:na], m[:soh][1,y+1,s] - 2 * storages[a].nCycle * m[:E_state][y] >= -storages[a].bounds.ub * storages[a].α_soc_max * (1-m[:invest_B][y])
        [y in 1:(ny-1), s in 1:ns, a in 1:na], m[:soh][1,y+1,s] - m[:soh][end,y,s] <= storages[a].bounds.ub * 2 * storages[a].nCycle
        [y in 1:(ny-1), s in 1:ns, a in 1:na], m[:soh][1,y+1,s] - m[:soh][end,y,s] >= -storages[a].bounds.ub * 2 * storages[a].nCycle

        #TODO verifier les Big M de mes 8 contraintes
        end)
    end
end
function add_technical_constraints!(m::Model, converters::Vector{AbstractConverter}, nh::Int64, ny::Int64, ns::Int64)
    if !isempty(converters)
        na = length(converters)
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, y in 1:ny, s in 1:ns, a in 1:na], m[:p_c][h,y,s,a]  <= m[:r_c][a]
        end)
    end
end
function add_technical_constraints!(m::Model, grids::Vector{AbstractGrid}, nh::Int64, ny::Int64, ns::Int64)
    if !isempty(grids)
        na = length(grids)
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, y in 1:ny, s in 1:ns, a in 1:na], m[:p_in][h,y,s,a]  <= grids[a].powerMax
        [h in 1:nh, y in 1:ny, s in 1:ns, a in 1:na], m[:p_out][h,y,s,a] <= grids[a].powerMax
        end)
    end
end



# Periodicity constraint
function add_periodicity_constraints!(m::Model, storages::Vector{AbstractStorage}, ns::Int64)
    # Storages
    if !isempty(storages)
        na = length(storages)
        @constraints(m, begin
        # Final states
        [s in 1:ns, a in 1:na], m[:soc][end,s,a]  >= m[:soc][1,s,a]
        end)
    end
end





# Power balance
function add_power_balance!(m::Model, mg::Microgrid, ω::Scenarios, type::DataType, nh::Int64, ns::Int64; ispnet::Bool=false)
    # !!! All the decision variables are defined positive !!!
    balance = AffExpr.(zeros(nh,ns))
    # Demands and generation
    if !ispnet
        for (k,a) in enumerate(mg.demands)
            if a.carrier isa type
                add_to_expression!.(balance, ω.demands[k].power[:,1,:])
            end
        end
        # Generation
        for (k,a) in enumerate(mg.generations)
            if a.carrier isa type
                add_to_expression!.(balance, .- m[:r_g][k] .* ω.generations[k].power[:,1,:])
            end
        end
    else
        for (k,a) in enumerate(mg.demands)
            if a.carrier isa type
                add_to_expression!.(balance, m[:p_d][:,:,k])
            end
        end
        # Generation
        for (k,a) in enumerate(mg.generations)
            if a.carrier isa type
                add_to_expression!.(balance, .- m[:p_g][:,:,k])
            end
        end
    end
    # Storages
    for (k,a) in enumerate(mg.storages)
        if a.carrier isa type
            add_to_expression!.(balance, m[:p_ch][:,:,k] .- m[:p_dch][:,:,k])
        end
    end
    # Converters
    for (k,a) in enumerate(mg.converters)
        if type == Electricity
            if a isa Heater
                add_to_expression!.(balance, m[:p_c][:,:,k])
            elseif a isa Electrolyzer
                add_to_expression!.(balance, m[:p_c][:,:,k])
            elseif a isa FuelCell
                add_to_expression!.(balance, .- m[:p_c][:,:,k])
            end
        elseif type == Heat
            if a isa Heater
                add_to_expression!.(balance, .- m[:p_c][:,:,k] * a.η_E_H)
            elseif a isa Electrolyzer
                add_to_expression!.(balance, .- m[:p_c][:,:,k] * a.η_E_H)
            elseif a isa FuelCell
                add_to_expression!.(balance, .- m[:p_c][:,:,k] * a.η_H2_H / a.η_H2_E)
            end
        elseif type == Hydrogen
            if a isa Electrolyzer
                add_to_expression!.(balance, .- m[:p_c][:,:,k] * a.η_E_H2)
            elseif a isa FuelCell
                add_to_expression!.(balance, m[:p_c][:,:,k] / a.η_H2_E)
            end
        end
    end
    # Grids
    for (k,a) in enumerate(mg.grids)
        if a.carrier isa type
            add_to_expression!.(balance, .- m[:p_in][:,:,k] + m[:p_out][:,:,k])
        end
    end
    # Energy balance constraint
    if type == Electricity
        @constraint(m, electricity, balance .<= 0.)
    elseif type == Heat
        @constraint(m, heat, balance .<= 0.)
    elseif type == Hydrogen
        @constraint(m, hydrogen, balance .== 0.)
    end
end



# Power balance multi year case
function add_power_balance!(m::Model, mg::Microgrid, ω::Scenarios, type::DataType, nh::Int64, ny::Int64, ns::Int64; ispnet::Bool=false)
    # !!! All the decision variables are defined positive !!!
    balance = AffExpr.(zeros(nh,ny,ns))
    # Demands and generation

    for y in 1:ny
        if !ispnet
            for (k,a) in enumerate(mg.demands)
                if a.carrier isa type
                    add_to_expression!.(balance, ω.demands[k].power[:,y,:])
                end
            end
            # Generation
            for (k,a) in enumerate(mg.generations)
                if a.carrier isa type
                    add_to_expression!.(balance, .- m[:r_g][k] .* ω.generations[k].power[:,1,:])
                end
            end
        else
            for (k,a) in enumerate(mg.demands)
                if a.carrier isa type
                    add_to_expression!.(balance, m[:p_d][:,y,:,k])
                end
            end
            # Generation
            for (k,a) in enumerate(mg.generations)
                if a.carrier isa type
                    add_to_expression!.(balance, .- m[:p_g][:,y,:,k])
                end
            end
        end


        # Storages
        for (k,a) in enumerate(mg.storages)
            if a.carrier isa type
                add_to_expression!.(balance, m[:p_ch][:,y,:,k] .- m[:p_dch][:,y,:,k])
            end
        end
        # Converters
        for (k,a) in enumerate(mg.converters)
            if type == Electricity
                if a isa Heater
                    add_to_expression!.(balance, m[:p_c][:,y,:,k])
                elseif a isa Electrolyzer
                    add_to_expression!.(balance, m[:p_c][:,y,:,k])
                elseif a isa FuelCell
                    add_to_expression!.(balance, .- m[:p_c][:,y,:,k])
                end
            elseif type == Heat
                if a isa Heater
                    add_to_expression!.(balance, .- m[:p_c][:,y,:,k] * a.η_E_H)
                elseif a isa Electrolyzer
                    add_to_expression!.(balance, .- m[:p_c][:,y,:,k] * a.η_E_H)
                elseif a isa FuelCell
                    add_to_expression!.(balance, .- m[:p_c][:,y,:,k] * a.η_H2_H / a.η_H2_E)
                end
            elseif type == Hydrogen
                if a isa Electrolyzer
                    add_to_expression!.(balance, .- m[:p_c][:,y,:,k] * a.η_E_H2)
                elseif a isa FuelCell
                    add_to_expression!.(balance, m[:p_c][:,y,:,k] / a.η_H2_E)
                end
            end
        end
        # Grids
        for (k,a) in enumerate(mg.grids)
            if a.carrier isa type
                add_to_expression!.(balance, .- m[:p_in][:,y,:,k] + m[:p_out][:,y,:,k])
            end
        end
    end
    # Energy balance constraint
    if type == Electricity
        @constraint(m, electricity, balance .<= 0.)
    elseif type == Heat
        @constraint(m, heat, balance .<= 0.)
    elseif type == Hydrogen
        @constraint(m, hydrogen, balance .== 0.)
    end
end


# Renewable share
function add_renewable_share!(m::Model, mg::Microgrid, ω::Scenarios, probabilities::Vector{Float64}, risk::AbstractRiskMeasure, nh::Int64, ns::Int64)
    total = zeros(ns)
    for (k,a) in enumerate(mg.demands)
        if a.carrier isa Electricity
            total .= total .+ sum(ω.demands[k].power[h,1,:] for h in 1:nh)
        elseif a.carrier isa Heat
            total .= total .+ sum(ω.demands[k].power[h,1,:] for h in 1:nh) ./ mg.converters[isin(mg.converters, Heater)[2]].η_E_H
        end
    end
    for (k,a) in enumerate(mg.grids)
        if a.carrier isa Electricity
            @expression(m, share[s in 1:ns], sum(m[:p_in][h,s,k] for h in 1:nh) - (1. - mg.parameters.renewable_share) * total[s])
        end
    end
    # Constraint according to CVaR
    @variables(m, begin
    ζ_s
    α_s[1:ns] >= 0.
    end)
    @constraints(m, begin
    [s in 1:ns], α_s[s] >= m[:share][s] - ζ_s
    ζ_s + 1 / (1 - beta(risk)) * sum(probabilities[s] * α_s[s] for s in 1:ns) <= 0.
    end)
end



# Renewable share in multi yeat context
function add_renewable_share!(m::Model, mg::Microgrid, ω::Scenarios, probabilities::Vector{Float64}, risk::AbstractRiskMeasure, nh::Int64, ny::Int64, ns::Int64)
    total = zeros(ns)
    for (k,a) in enumerate(mg.demands)
        if a.carrier isa Electricity
            total .= total .+ sum(ω.demands[k].power[h,y,:] for h in 1:nh for y in 1:ny)
        elseif a.carrier isa Heat
            total .= total .+ sum(ω.demands[k].power[h,y,:] for h in 1:nh for y in 1:ny) ./ mg.converters[isin(mg.converters, Heater)[2]].η_E_H
        end
    end
    for (k,a) in enumerate(mg.grids)
        if a.carrier isa Electricity
            @expression(m, share[s in 1:ns], sum(m[:p_in][h,y,s,k] for h in 1:nh for y in 1:ny) - (1. - mg.parameters.renewable_share) * total[s])
        end
    end
    # Constraint according to CVaR
    @variables(m, begin
    ζ_s
    α_s[1:ns] >= 0.
    end)
    @constraints(m, begin
    [s in 1:ns], α_s[s] >= m[:share][s] - ζ_s
    ζ_s + 1 / (1 - beta(risk)) * sum(probabilities[s] * α_s[s] for s in 1:ns) <= 0.
    end)
end




# Objective
function add_design_objective!(m::Model, mg::Microgrid, ω::Scenarios, probabilities::Vector{Float64}, risk::AbstractRiskMeasure, nh::Int64, ns::Int64)
    # CAPEX
    capex = compute_capex(m, mg, ω)
    # OPEX
    opex = compute_opex(m, mg, ω, nh, ns)
    # Objective according to the CVaR
    @variables(m, begin
    ζ_o
    α_o[1:ns] >= 0.
    end)
    @constraint(m, [s in 1:ns], α_o[s] >= capex + opex[s] - ζ_o)
    @objective(m, Min, ζ_o + 1 / (1 - beta(risk)) * sum(probabilities[s] * α_o[s] for s in 1:ns))
end


# Capex
function compute_capex(m::Model, mg::Microgrid, ω::Scenarios)
    cost = AffExpr(0.)
    # Generations
    for (k,a) in enumerate(mg.generations)
        add_to_expression!(cost, Γ(mg.parameters.τ, a.lifetime) * ω.generations[k].cost[1] * m[:r_g][k])
    end
    # Storages
    for (k,a) in enumerate(mg.storages)
        add_to_expression!(cost, Γ(mg.parameters.τ, a.lifetime) * ω.storages[k].cost[1] * m[:r_sto][k])
    end
    # Converters
    for (k,a) in enumerate(mg.converters)
        add_to_expression!(cost, Γ(mg.parameters.τ, a.lifetime) * ω.converters[k].cost[1] * m[:r_c][k])
    end
    return cost
end

# Grids
function compute_opex(m::Model, mg::Microgrid, ω::Scenarios, nh::Int64, ns::Int64)
    cost = AffExpr.(zeros(ns))
    for (k,a) in enumerate(mg.grids)
        add_to_expression!.(cost, sum((m[:p_in][h,:,k] .* ω.grids[k].cost_in[h,1,:] .- m[:p_out][h,:] .* ω.grids[k].cost_out[h,1,:]) .* mg.parameters.Δh  for h in 1:nh))
    end
    return cost
end


#Multiple year dynamic case
# Objective
function add_design_objective!(m::Model, mg::Microgrid, ω::Scenarios, probabilities::Vector{Float64}, risk::AbstractRiskMeasure, nh::Int64, ny::Int64, ns::Int64)
    # CAPEX
    capex = compute_capex(m, mg, ω, ny)
    # OPEX
    opex = compute_opex(m, mg, ω, nh, ny, ns)

    salvage = compute_salvage(m, mg, ω, ny, ns)
    # Objective according to the CVaR


    @objective(m, Min,  capex + opex[1] - salvage[1])
end


# Capex
function compute_capex(m::Model, mg::Microgrid, ω::Scenarios, ny::Int64)
    cost = AffExpr(0.)
    for y in 1:ny
        # Generations
        for (k,a) in enumerate(mg.generations)
            add_to_expression!(cost, ω.generations[k].cost[y] * m[:r_g][k])
        end
        # Storages
        for (k,a) in enumerate(mg.storages)
            add_to_expression!(cost, ω.storages[k].cost[y] * m[:r_sto][y,k])
        end
        # Converters
        # for (k,a) in enumerate(mg.converters)
        #     add_to_expression!(cost, ω.converters[k].cost[y] * m[:r_c][y,k])
        # end
    end
    return cost
end

# Grids
function compute_opex(m::Model, mg::Microgrid, ω::Scenarios, nh::Int64, ny::Int64, ns::Int64)
    cost = AffExpr.(zeros(ns))
    for (k,a) in enumerate(mg.grids)
        add_to_expression!.(cost, sum((m[:p_in][h,y,:,k] .* ω.grids[k].cost_in[h,y,:] .- m[:p_out][h,y,:] .* ω.grids[k].cost_out[h,y,:]) .* mg.parameters.Δh  for h in 1:nh for y in 1:ny))
    end
    return cost
end

function compute_salvage(m::Model, mg::Microgrid, ω::Scenarios, ny::Int64, ns::Int64)
    salvage =  AffExpr.(zeros(ns))
    K = ω.storages[1].cost[ny] * m[:E_state][end]
    add_to_expression!.(salvage, K * m[:soh][end,end,:])

    return salvage
end
