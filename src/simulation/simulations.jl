#=
    This file includes all the funtions needed to simulate the DES
 =#

 mutable struct Options
     mode::String
     firstyear::Bool #simulate
     Options(; mode = "multithreads", firstyear = false) = new(mode, firstyear)
 end

# Simulate function
function simulate!(mg::Microgrid,
                   controller::AbstractController,
                   designer::AbstractDesigner,
                   ω_simu::Scenarios;
                   options::Options = Options())


    #TODO Maybe need to reset some parameters like calendar aging on electochimical batteries
    # Parameters
    ns = mg.parameters.ns

    if options.mode == "serial"
        # We simulate over the horizon for all the scenarios
        @showprogress for s in 1:ns
            simulate!(s, mg, controller, designer, ω_simu, options)
        end

    elseif options.mode == "multicores"
        # Init
        s = 0
        # We simulate over the horizon for all the scenarios in parallel
        @sync begin
            for p in workers()
                @async begin
                    while true
                        # Increment scenario
                        s += 1
                        # Break if ns is reached...
                        if s > ns
                            break
                        end
                        # Execute the function
                        remotecall_fetch(simulate!, p, s, mg, controller, designer, ω_simu, options)
                    end
                end
            end
        end
    elseif options.mode == "distributed"
        # We simulate over the horizon for all the scenarios in parallel using the distributed macro
        @sync @distributed for s in 1:ns
            simulate!(s, mg, controller, designer, ω_simu, options)
        end
    elseif options.mode == "multithreads"
        # We simulate over the horizon for all the scenarios in parallel using the distributed macro
        Threads.@threads for s in 1:ns
            simulate!(s, mg, controller, designer, ω_simu, options)
        end
    else
        println("Unknown mode... Please chose between 'serial', 'multicores', 'distributed' or 'multithreads'.")
    end
end
function simulate!(s::Int64,
                   mg::Microgrid,
                   controller::AbstractController,
                   designer::AbstractDesigner,
                   ω_simu::Scenarios,
                   options::Options)

    # Parameters
    nh = mg.parameters.nh
    ny = mg.parameters.ny



    # We simulate over the horizon for a single scenario
    for y in 1:ny
        simulate!(y, s, mg, controller, designer, ω_simu, options)
    end
end
function simulate!(y::Int64,
                   s::Int64,
                   mg::Microgrid,
                   controller::AbstractController,
                   designer::AbstractDesigner,
                   ω_simu::Scenarios,
                   options::Options)

    # Parameters
    nh = mg.parameters.nh

    if y == 1 && !options.firstyear
        # Update investment informations
        update_investment_informations!(y, s, mg, ω_simu)

        # Compute investment decision variables
        compute_investment_decisions!(y, s, mg, designer)

        # Compute investment dynamics
        compute_investment_dynamics!(y, s, mg, designer)

        #update grid prices
        update_grid_cost_informations!(y, s, mg, ω_simu)
    else
        for h in 1:nh
            simulate!(h, y, s, mg, controller, designer, ω_simu, options)
        end

        # Update investment informations
        update_investment_informations!(y, s, mg, ω_simu)

        # Compute investment decision variables
        compute_investment_decisions!(y, s, mg, designer)

        # Compute investment dynamics
        compute_investment_dynamics!(y, s, mg, designer)

        #update grid prices
        update_grid_cost_informations!(y, s, mg, ω_simu)


    end
end
function simulate!(h::Int64,
                   y::Int64,
                   s::Int64,
                   mg::Microgrid,
                   controller::AbstractController,
                   designer::AbstractDesigner,
                   ω_simu::Scenarios,
                   options::Options)

    # Update operation informations
    update_operation_informations!(h, y, s, mg, ω_simu)

    # Compute operation decision variables
    compute_operation_decisions!(h, y, s, mg, controller)

    # Compute operation dynamics for each converter and storage in mg
    compute_operation_dynamics!(h, y, s, mg, controller)

    # Power balance constraint checked for each node
    #Si désactivé rien n'est puisé sur la grid 
    compute_power_balances!(h, y, s, mg)
end
