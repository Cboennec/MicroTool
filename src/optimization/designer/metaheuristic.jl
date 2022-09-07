#=
    Designer based on a metaheuristic
=#

mutable struct MetaheuristicOptions
    method::Metaheuristics.AbstractMetaheuristic
    iterations::Int64
    multithreads::Bool
    controller::AbstractController
    isnpv::Bool
    reducer::AbstractScenariosReducer
    objective_risk::AbstractRiskMeasure
    share_risk::AbstractRiskMeasure
    lpsp_risk::AbstractRiskMeasure
    lpsp_tol::Float64
    reopt::Bool
    read_reduction::Union{String, Nothing}
    write_reduction::Union{String, Nothing}

    MetaheuristicOptions(; method = Metaheuristics.Clearing(),
                           iterations = 50,
                           multithreads = false,
                           controller = RBC(),
                           isnpv = false,
                           reducer = FeatureBasedReducer(),
                           objective_risk = Expectation(),
                           share_risk = Expectation(),
                           lpsp_risk = WorstCase(),
                           lpsp_tol = 1e-3,
                           reopt = false,
                           read_reduction = nothing,
                           write_reduction = nothing) =
                           new(method, iterations, multithreads, controller, isnpv, reducer, objective_risk, share_risk, lpsp_risk, lpsp_tol, reopt, read_reduction, write_reduction)

end

mutable struct Metaheuristic <: AbstractDesigner
    options::MetaheuristicOptions
    decisions::NamedTuple
    results::Metaheuristics.MetaheuristicResults
    history::AbstractScenarios

    Metaheuristic(; options = MetaheuristicOptions()) = new(options)
end

# Objective functions
function fobj(decisions::Array{Float64,1}, mg::Microgrid, designer::Metaheuristic, ω::Scenarios, probabilities::Array{Float64})
    # Paramters
    nh, ny, ns = size(ω.demands[1].power)
    λ1 = λ2 = λ3 = 1e6


    # Initialize mg
    mg_m = deepcopy(mg)
    mg_m.parameters = GlobalParameters(nh, ny, ns)



    # Initialize controller
    controller_m = initialize_controller!(mg_m, designer.options.controller, ω)

    # Initialize with the manual designer
    designer_m = initialize_designer!(mg_m, Manual(generations = [decisions[1:length(mg_m.generations)]...], storages = [decisions[length(mg_m.generations)+1:length(mg_m.generations)+length(mg_m.storages)]...], converters = [decisions[end-length(mg_m.converters)+1:end]...]), ω)

    # Simulate
    simulate!(mg_m, controller_m, designer_m, ω)



    # Metrics
    metrics = Metrics(mg_m, designer_m)


    # Share constraint
    #share = max(0., mg.parameters.renewable_share - conditional_value_at_risk([reshape(metrics.renewable_share[2:ny, 1:ns], :, 1)...],  probabilities,  designer.options.share_risk))
    share = max(0., mg.parameters.renewable_share - conditional_value_at_risk([reshape(metrics.renewable_share[1, 1:ns], :, 1)...],  probabilities,  designer.options.share_risk))

    # LPSP constraint for the heat
    metrics.lpsp.heat isa Nothing ? lpsp = 0. : lpsp = max(0., conditional_value_at_risk([reshape(metrics.lpsp.heat[2:ny, 1:ns], :, 1)...], probabilities,  designer.options.lpsp_risk) - designer.options.lpsp_tol)

    # SoC constraint for the seasonal storage
    soc_seasonal = 0.
    for a in mg_m.storages
        if a isa H2Tank
            soc_seasonal += sum(max(0., a.soc[1,y,s] - a.soc[end,y,s]) for y in 2:ny, s in 1:ns)
        end
    end


    # Objective - Algortihm find the maximum
    if designer.options.isnpv
        # NPV
        npv = conditional_value_at_risk([metrics.npv.total...], probabilities, designer.options.objective_risk)
        return npv - λ1 * share - λ2 * lpsp - λ3 * soc_seasonal
    else
        # Equivalent annual cost
        eac = conditional_value_at_risk([metrics.eac.total...], probabilities, designer.options.objective_risk)
        return - eac - λ1 * share - λ2 * lpsp - λ3 * soc_seasonal
    end
end

### Offline
function initialize_designer!(mg::Microgrid, designer::Metaheuristic, ω::Scenarios)
    # Preallocate and assigned values
    preallocate!(mg, designer)

    # Scenario reduction from the optimization scenario pool
    if designer.options.isnpv
        println("Starting scenario reduction...")
        ω_reduced, probabilities = reduce(designer.options.reducer, ω)
    else
        if isa(designer.options.read_reduction, Nothing)
            println("Starting scenario reduction...")
            ω_reduced, probabilities = reduce(designer.options.reducer, ω)
            # Saving
            if !isa(designer.options.write_reduction, Nothing)
                save(designer.options.write_reduction, "scenarios", ω_reduced, "probabilities", probabilities)
            end
        else
            println("Reading scenario reduction from file...")
            ω_reduced = load(designer.options.read_reduction, "scenarios")
            probabilities = load(designer.options.read_reduction, "probabilities")
        end
        # Repeat to simulate 2 years
        #ω_reduced = repeat(ω_reduced, 1, 4, 1)
    end

    # Bounds
    lb, ub = set_bounds(mg)

    println(probabilities)



    # Optimize
    designer.results = Metaheuristics.optimize(lb, ub,
                                               designer.options.method,
                                               options = Metaheuristics.Options(iterations = designer.options.iterations, multithreads = designer.options.multithreads)
    ) do decisions
        fobj(decisions, mg, designer, ω_reduced, probabilities)
      end

    # Assign values
    for k in 1:length(mg.generations)
        designer.decisions.generations[k][1,:] .= designer.results.minimizer[k]
    end
    for k in 1:length(mg.storages)
        designer.decisions.storages[k][1,:] .= designer.results.minimizer[length(mg.generations)+k]
    end
    for k in 1:length(mg.converters)
        designer.decisions.converters[k][1,:] .= designer.results.minimizer[end-length(mg.converters)+k]
    end

    # Save history for online optimization
    designer.history = ω_reduced

    return designer
end

### Online
# Loi de gestion d'investissement dynamique
# ex : remplacer la batterie par une equivalente à partir d'un seuil défini
# regarder dans le papier sur investissement dynamique
function compute_investment_decisions!(y::Int64, s::Int64, mg::Microgrid, designer::Metaheuristic)

    if mg.storages[1].soh[end,y,s] <= mg.storages[1].SoH_threshold
        designer.decisions.storages[1][y,s] =  designer.storages[1]
    end
    return nothing
end

### Utils
function set_bounds(mg::Microgrid)
    # Initialization
    lb, ub = [], []
    # Generations
    for a in mg.generations
        push!(lb, a.bounds.lb)
        push!(ub, a.bounds.ub)
    end
    # Storages
    for a in mg.storages
        push!(lb, a.bounds.lb)
        push!(ub, a.bounds.ub)
    end
    # Converters
    for a in mg.converters
        push!(lb, a.bounds.lb)
        push!(ub, a.bounds.ub)
    end
    return lb, ub
end
