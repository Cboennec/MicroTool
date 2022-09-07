#=
    Rule based controller
=#
mutable struct RBCOptions
    policy_selection::Int64

    RBCOptions(; policy_selection = 1) = new(policy_selection)
end

mutable struct RBC <: AbstractController
    options::RBCOptions
    decisions::NamedTuple
    history::AbstractScenarios

    RBC(; options = RBCOptions()) = new(options)
end

### Policies
function π_1(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::RBC)
    # Utils to simplify the writting
    Δh = mg.parameters.Δh
    liion, tes, h2tank = mg.storages[1], mg.storages[2], mg.storages[3]
    heater, elyz, fc = mg.converters[1], mg.converters[2], mg.converters[3]

    # Net power elec
    p_net_E = mg.demands[1].carrier.power[h,y,s] - mg.generations[1].carrier.power[h,y,s]

    # Liion
    _, _, u_liion = compute_operation_dynamics(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]), p_net_E, Δh)

    if p_net_E < 0.
        # Elyz
        _, u_elyz_E, elyz_H, elyz_H2 = compute_operation_dynamics(elyz, (powerMax = elyz.powerMax[y,s], soh = elyz.soh[h,y,s]), p_net_E - u_liion, Δh)
        # H2 tank
        _, u_h2tank = compute_operation_dynamics(h2tank, (Erated = h2tank.Erated[y,s], soc = h2tank.soc[h,y,s]), -elyz_H2, Δh)
        # Test H2
        elyz_H2 == - u_h2tank ? nothing : u_elyz_E = elyz_H = elyz_H2 = u_h2tank = 0.
        # FC
        u_fc_E, fc_H, fc_H2 = 0., 0., 0.
        # Heater
        u_heater_E, heater_H = compute_operation_dynamics(heater, (powerMax = heater.powerMax[y,s],), p_net_E - u_liion - u_elyz_E, Δh)
    else
        # FC
        _, u_fc_E, fc_H, fc_H2 = compute_operation_dynamics(fc, (powerMax = fc.powerMax[y,s], soh = fc.soh[h,y,s]), p_net_E - u_liion, Δh)
        # H2 tank
        _, u_h2tank = compute_operation_dynamics(h2tank, (Erated = h2tank.Erated[y,s], soc = h2tank.soc[h,y,s]), -fc_H2, Δh)
        # Test H2
        fc_H2 == - u_h2tank ? nothing : u_fc_E = fc_H = fc_H2 = u_h2tank = 0.
        # Elyz
        u_elyz_E, elyz_H, elyz_H2 = 0., 0., 0.
        # Heater
        u_heater_E, heater_H = 0., 0.
    end

    # Net heating power post H2
    p_net_H = mg.demands[2].carrier.power[h,y,s] - fc_H - elyz_H - heater_H

    # TES
    _, u_tes = compute_operation_dynamics(tes, (Erated = tes.Erated[y,s], soc = tes.soc[h,y,s]), p_net_H, Δh)

    # Heater
    if p_net_H < 0.
        _u_heater_E = 0.
    else
        _u_heater_E, _ = compute_operation_dynamics(heater, (powerMax = heater.powerMax[y,s],), - (p_net_H - u_tes) / heater.η_E_H, Δh)
    end

    # Store values
    controller.decisions.storages[1][h,y,s] = u_liion
    controller.decisions.storages[2][h,y,s] = u_tes
    controller.decisions.storages[3][h,y,s] = u_h2tank
    controller.decisions.converters[1][h,y,s]  = u_heater_E + _u_heater_E
    controller.decisions.converters[2][h,y,s] = u_elyz_E
    controller.decisions.converters[3][h,y,s] = u_fc_E
end
function π_2(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::RBC)
    controller.decisions.storages[1][h,y,s] = mg.demands[1].carrier.power[h,y,s] - mg.generations[1].carrier.power[h,y,s]
end
function π_3(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::RBC)
    # Compute the heater electrical power based on the simple model
    controller.decisions.converters[1][h,y,s] = - max(min(mg.demands[2].carrier.power[h,y,s] / mg.converters[1].η_E_H, mg.converters[1].powerMax[y,s]), 0.)
    # Compute the liion decision from the power balance
    controller.decisions.storages[1][h,y,s] = mg.demands[1].carrier.power[h,y,s] - mg.generations[1].carrier.power[h,y,s] - controller.decisions.converters[1][h,y,s]
end





### Offline
function initialize_controller!(mg::Microgrid, controller::RBC, ω::AbstractScenarios)
    # Preallocation
    preallocate!(mg, controller)

    return controller
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::RBC)
    # Chose policy
    if controller.options.policy_selection == 1
        return π_1(h, y, s, mg, controller)
    elseif controller.options.policy_selection == 2
        return π_2(h, y, s, mg, controller)
    elseif controller.options.policy_selection == 3
        return π_3(h, y, s, mg, controller)
    else
        println("Policy not defined !")
    end
end
