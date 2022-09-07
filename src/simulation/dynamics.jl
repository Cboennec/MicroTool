#=
    This file includes all the funtions needed to compute the operation
    and investment dynamics
 =#

function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::AbstractController)
    # Converters
    for (k, a) in enumerate(mg.converters)
        compute_operation_dynamics!(h, y, s, a, controller.decisions.converters[k][h,y,s], mg.parameters.Δh)
    end
    # Storage
    for (k, a) in enumerate(mg.storages)
        compute_operation_dynamics!(h, y, s, a, controller.decisions.storages[k][h,y,s], mg.parameters.Δh)
    end
end
function compute_investment_dynamics!(y::Int64, s::Int64, mg::Microgrid, designer::AbstractDesigner)
     # Generations
     for (k, a) in enumerate(mg.generations)
        compute_investment_dynamics!(y, s, a, designer.decisions.generations[k][y,s])
     end
     # Storages
     for (k, a) in enumerate(mg.storages)
         compute_investment_dynamics!(y, s, a, designer.decisions.storages[k][y,s])
     end
     # Converters
     for (k, a) in enumerate(mg.converters)
         compute_investment_dynamics!(y, s, a, designer.decisions.converters[k][y,s])
     end
end
