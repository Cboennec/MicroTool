#=
    Model predictive control controller
=#

mutable struct OLFCOptions
    solver
    reducer::AbstractScenariosReducer
    generator::AbstractScenariosGenerator
    horizon::Int64
    nscenarios::Int64
    seasonal_targets::Union{Array{Float64,2}, Nothing}
    final_cost::Vector{Float64}

    OLFCOptions(; solver = Cbc,
                  reducer = FeatureBasedReducer(),
                  generator = MarkovGenerator(),
                  horizon = 24,
                  nscenarios = 1,
                  seasonal_targets = nothing,
                  final_cost = [1e-3, 1e-3, 1e-2]) = new(solver, reducer, generator, horizon, nscenarios, seasonal_targets, final_cost)
end

mutable struct OLFC <: AbstractController
    options::OLFCOptions
    generations::Vector{Float64}
    storages::Vector{Float64}
    converters::Vector{Float64}
    decisions::NamedTuple
    model::JuMP.Model
    history::AbstractScenarios
    OLFC(; options = OLFCOptions(),
           generations = [0.],
           storages = [0.],
           converters = []) =
           new(options, generations, storages, converters)
end

### Model
function build_model(mg::Microgrid, controller::OLFC, ω::Scenarios)
    # Sets
    nh, ns = controller.options.horizon, controller.options.nscenarios
    # Initialize
    m = Model(controller.options.solver.Optimizer)
    set_optimizer_attribute(m,"CPX_PARAM_SCRIND", 0)
    # Add investment variables
    add_investment_decisions!(m, mg.generations)
    add_investment_decisions!(m, mg.storages)

    add_investment_decisions!(m, mg.converters)
    # Fix their values
#    fix_investment_decisions!(m, controller.generations, controller.storages, controller.converters)
    # Add operation decision variable with the non-anticipative structure
    add_operation_decisions!(m, mg.storages, nh, 1)
    if !isempty(controller.converters)
        add_operation_decisions!(m, mg.converters, nh, 1)
    end
    # Add demand and generation variables to be fixed online, along with recourse variable (grid)
    add_operation_decisions!(m, mg.generations, nh, ns)
    add_operation_decisions!(m, mg.demands, nh, ns)
    add_operation_decisions!(m, mg.grids, nh, ns)
    # Add technical constraints
    add_technical_constraints!(m, mg.storages, mg.parameters.Δh, nh, 1)
    if !isempty(controller.converters)
#        add_technical_constraints!(m, mg.converters, nh, 1)
    end
#    add_technical_constraints!(m, mg.grids, nh, ns)
    # Add power balance constraints
    add_power_balance!(m, mg, ω, Electricity, nh, ns, ispnet = true)
#    add_power_balance!(m, mg, ω, Heat, nh, ns, ispnet = true)
#    add_power_balance!(m, mg, ω, EnergyCarrier, nh, ns, ispnet = true)
    return m
end

### Offline
function initialize_controller!(mg::Microgrid, controller::OLFC, ω::Scenarios)
    # Preallocate
    preallocate!(mg, controller)
    # Build model
    println("Building the model...")
    controller.model = build_model(mg, controller, ω)
    # Scenario reduction
    println("Starting scenario reduction...")
    ω_reduced = reduce(controller.options.reducer, ω)[1]
    # Compute markov chain for scenario generation
    println("Initializing scenario generator...")
    controller.options.generator = initialize_generator!(controller.options.generator, [a for a in ω_reduced.generations]..., [a for a in ω_reduced.demands]...)
    # History
    controller.history = ω
    return controller
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::OLFC)
    # Forecast window
    window = h:min(mg.parameters.nh, h + controller.options.horizon - 1)
    # Compute forecast
    demands, generations, grids, probabilities = compute_forecast(h, y, s, mg, controller, window)
    # Fix model variables
    fix_operation_variables!(controller, demands, generations, [a.soc[h,y,s] * a.Erated[y,s] for a in mg.storages])
    # Add objective
    add_objective!(controller, mg, grids, probabilities, window)
    # Optimize
    optimize!(controller.model)
    # Assign controller values
    for k in 1:length(mg.storages)
        controller.decisions.storages[k][h,y,s] = value(controller.model[:p_dch][1,1,k] - controller.model[:p_ch][1,1,k])
    end
    for (k,a) in enumerate(mg.converters)
        if a isa Heater
            controller.decisions.converters[k][h,y,s] = - value(controller.model[:p_c][1,1,k])
        elseif a isa Electrolyzer
            controller.decisions.converters[k][h,y,s] = - value(controller.model[:p_c][1,1,k])
        elseif a isa FuelCell
            controller.decisions.converters[k][h,y,s] = value(controller.model[:p_c][1,1,k])
        end
    end
end
# OLFC forecast
function compute_forecast(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::OLFC, window::UnitRange{Int64})
    # To make the writting easier...
    nh, ns = controller.options.horizon, controller.options.nscenarios
    # Initial state
    s0 = vcat([a.carrier.power[h,y,s] / a.powerMax[y,s] for a in mg.generations], [a.carrier.power[h,y,s] for a in mg.demands])
    # Initial timestamp
    t0 = mg.demands[1].timestamp[h,y,s]
    # Demand and production forecast using the scenario generator
    forecast, probabilities = generate(controller.options.generator, s0, t0, nh, ny = ns)
    # Unpack the forecast vector
    generations, demands = [a.powerMax[y,s] * forecast[k] for (k,a) in enumerate(mg.generations)], forecast[length(mg.generations)+1:end]
    # Operating cost - unchanged along the horizon. We add zeros to guarantee a constant vector size
    grids = [(cost_in = vcat(a.cost_in[window,y,s], zeros(nh - length(window))), cost_out = vcat(a.cost_out[window,y,s], zeros(nh - length(window)))) for a in mg.grids]

    return demands, generations, grids, probabilities
end
# Fix OLFC variables
function fix_operation_variables!(controller::OLFC, demands::Vector{Array{Float64,3}}, generations::Vector{Array{Float64,3}}, soc_ini::Vector{Float64})
    # Fix demand values
    for (k,a) in enumerate(demands)
        fix.(controller.model[:p_d][:,:,k], a[:,:,1])
    end
    # Fix generation values
    for (k,a) in enumerate(generations)
        fix.(controller.model[:p_g][:,:,k], a[:,:,1])
    end
    # Fix initial SoC values
    for (k,a) in enumerate(soc_ini)
        set_normalized_coefficient.(controller.model[:soc_ini][:,k], controller.model[:r_sto][k], 0.)
        set_normalized_rhs.(controller.model[:soc_ini][:,k], a)
    end
end
# Add objective
function add_objective!(controller::OLFC, mg::Microgrid, tariff::Vector{NamedTuple{(:cost_in, :cost_out), Tuple{Vector{Float64},Vector{Float64}}}}, probabilities::Array{Float64,2}, window::UnitRange{Int64})
    # To make the writting easier...
    nh, ns = controller.options.horizon, controller.options.nscenarios
    # Initialize
    cost = AffExpr(0.)
    # Opex
    for (k,a) in enumerate(mg.grids)
        add_to_expression!(cost, sum(probabilities[s] * sum((controller.model[:p_in][h,s,k] * tariff[k].cost_in[h] - controller.model[:p_out][h,s,k] * tariff[k].cost_out[h]) * mg.parameters.Δh for h in 1:nh) for s in 1:ns))
    end
    # Add final cost
    for (k,a) in enumerate(mg.storages)
        if controller.options.seasonal_targets isa Nothing
            add_to_expression!(cost, - controller.options.final_cost[k] * controller.model[:soc][end,1,k])
        else
            if a isa H2Tank
                add_to_expression!(cost, - controller.options.final_cost[k] * (controller.model[:soc][end,1,k] - controller.options.seasonal_targets[window]))
            else
                add_to_expression!(cost, - controller.options.final_cost[k] * controller.model[:soc][end,1,k])
            end
        end
    end
    @objective(controller.model, Min, cost)
end
