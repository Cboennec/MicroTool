# This file includes all the plot functions for the MG

using Pandas

function plot_operation(mg::Microgrid ; y=2, s=1, smooth = false, xdisplay = "hours", ci = false)
    # Seaborn configuration
    Seaborn.set_theme(context="notebook", style="ticks", palette="muted", font="serif", font_scale=1.5)


    # Parameters
    nh = mg.parameters.nh
    Δh = mg.parameters.Δh
    hours = range(1, length = nh * length(y), step = Δh) / Δh



    #If there is no converter let's say we have only one type of energy carrier (can theoriticaly be false)
    if length(mg.converters) == 0
        energy_carriers = enumerate([typeof(mg.generations[1].carrier)])
    else #Else we enumerate what type of carrier we have in the converters
        energy_carriers_list = []
        for conv in mg.converters
            for carrier in conv.carrier
                push!(energy_carriers_list, carrier)
            end
        end

        energy_carriers = enumerate([unique((typeof(a) for a in energy_carriers_list))])
    end

    # Plots
    # Powers
    f = figure("Powers")
    f.subplotpars.hspace = 0.32
    for (i, type) in energy_carriers
        i == 1 ? subplot(3, 1, i, title = string(type)) : subplot(3, 1, i, sharex = f.axes[1], title = string(type))
        # Demands
        for (k, a) in enumerate(mg.demands)
            if a.carrier isa type

                Seaborn.plot(hours, vec(a.carrier.power[:,y,s]), label = string("Demand : ",  typeof(a.carrier)))
            end
        end
        # Generations
        for (k, a) in enumerate(mg.generations)
            if a.carrier isa type
                Seaborn.plot(hours, vec(a.carrier.power[:,y,s]), label = string("Generation : ", typeof(a)))
            end
        end
        # Storages
        for (k, a) in enumerate(mg.storages)
            if a.carrier isa type
                Seaborn.plot(hours, vec(a.carrier.power[:,y,s]), label = string("Storage : ", typeof(a)))
            end
        end
        # Converters
        for (k, a) in enumerate(mg.converters)
            for c in a.carrier
                if c isa type
                    Seaborn.plot(hours, vec(c.power[:,y,s]), label = string("Converter : ", typeof(a)))
                end
            end
        end
        for (k, a) in enumerate(mg.grids)
            if a.carrier isa type
                Seaborn.plot(hours,  vec(a.carrier.power[:,y,s]), label = string("Grids : ", typeof(a)))
            end
        end
        legend()
    end
    # State of charge
    figure("State-of-charge")
    for (k, a) in enumerate(mg.storages)
        k == 1 ? subplot(length(mg.storages), 1, k) : subplot(length(mg.storages), 1, k, sharex = f.axes[1])
        Seaborn.plot(hours, vec(a.soc[1:end-1, y, s]), label = string("Storage : ", typeof(a)))
        legend()
    end


    #State of health

    x_val = []
    y_val = []
    x_lab = ""
    soc = ""
    x_ticks_lab = []
    x_ticks = []

    for s in 1:mg.parameters.ns
        figure("State-of-health")

        mg.storages[1].soc_model == "linear" ? soc = "lin" : (mg.storages[1].soc_model == "vermeer" ? soc = "ver" : soc = "t-d")

        if xdisplay == "years"
            y_values = vec(mg.storages[1].soh[:, y[1]:mg.parameters.ny, s])[1:24:((y[end] - y[1] +1) * (mg.parameters.nh+1) )]
            x_values = 1:((nh / 24) * (y[end] - y[1] +1) +1)
            x_lab = "Years"
            x_ticks_lab = [string(n) for n in 0:y[end-1]]
            x_ticks = [365 * n for n in 0:y[end-1]]

        elseif xdisplay == "hours"
            y_values = vec(mg.storages[1].soh[:, y[1]:mg.parameters.ny, s])[1:((y[end] - y[1] +1) * (mg.parameters.nh+1) )]
            x_values = 1:((nh+1) * (y[end] - y[1] +1))
            x_lab = "Hours"
        elseif xdisplay == "days"
            y_values = vec(mg.storages[1].soh[:, y[1]:mg.parameters.ny, s])[1:24:((y[end] - y[1] +1) * (mg.parameters.nh+1) )]
            x_values = 1:((nh / 24) * (y[end] - y[1] +1) +1)
            x_lab = "Days"

        end

        if ci
            x_val = vcat(x_val, x_values)
            y_val = vcat(y_val, y_values)
        else


            if (mg.storages[1] isa Liion_electro_chimique || mg.storages[1] isa Liion_rainflow) && smooth
                plot_soh = Seaborn.plot(1:(mg.parameters.nh/mg.storages[1].update_by_year):((y[end] - y[1] +1) * (mg.parameters.nh+1) ) ,
                vec(mg.storages[1].soh[:, y[1]:mg.parameters.ny, s])[1:(convert(Int64,mg.parameters.nh/mg.storages[1].update_by_year)):((y[end] - y[1] +1) * (mg.parameters.nh+1) )],
                label = string("soh : ", typeof(mg.storages[1]), " ", mg.storages[1].couplage, " ", soc ))
            else
                plot_soh = Seaborn.plot(x_values ,
                y_values,
                label = string("soh : ", typeof(mg.storages[1]), " ", mg.storages[1].couplage, " ", soc, ", s",s ))
            end

            xlabel(x_lab, fontsize = 20)
            ylabel("SoH", fontsize = 20)
            legend()
        end
    end

    if ci
        df = Pandas.DataFrame(Dict(:x_val=>x_val, :y_val=>y_val))
        pl = Seaborn.lineplot(data=df, x="x_val", y="y_val", label = string("soh : ", typeof(mg.storages[1]), " ", mg.storages[1].couplage, " ", soc ), ci = 95)
    	xlabel(x_lab, fontsize = 20)
        ylabel("SoH", fontsize = 20)
        pl.set_xticks(x_ticks)
        pl.set_xticklabels(x_ticks_lab)


        legend()

    end


    #voltage
    if mg.storages[1].soc_model == "tremblay_dessaint" 
        figure("voltage")
        Seaborn.plot(1:((mg.parameters.ny -1) * (mg.parameters.nh+1) ), vec(mg.storages[1].voltage[:, 2:mg.parameters.ny, s]), label = string("voltage : ", typeof(mg.storages[1])))
        legend()
    end





end
function plot_operation(mg::Microgrid, controller::AbstractController; y=2, s=1)
    # Seaborn configuration
    Seaborn.set_theme(context="notebook", style="ticks", palette="muted", font="serif", font_scale=1.5)

    # Parameters
    nh = mg.parameters.nh
    Δh = mg.parameters.Δh
    hours = range(1, length = nh, step = Δh) / Δh

    # Plots
    # Powers
    f = figure("Powers")
    for (i, type) in enumerate([typeof(Electricity()), typeof(Heat()), typeof(Hydrogen())])
        i == 1 ? subplot(3, 1, i) : subplot(3, 1, i, sharex = f.axes[1])
        # Demands
        for (k, a) in enumerate(mg.demands)
            if a.carrier isa type
                plot(hours, a.carrier.power[:,y,s], label = string("Demand ", k))
            end
        end
        # Generations
        for (k, a) in enumerate(mg.generations)
            if a.carrier isa type
                plot(hours, a.carrier.power[:,y,s], label = string("Generation ", k))
            end
        end
        # Storages
        for (k, a) in enumerate(mg.storages)
            if a.carrier isa type
                plot(hours, controller.decisions.storages[k][:,y,s], label = string("Storage ", k))
            end
        end
        # Converters
        for (k, a) in enumerate(mg.converters)
            for c in a.carrier
                if c isa type
                    plot(hours, controller.decisions.converters[k][:,y,s], label = string("Converter ", k))
                end
            end
        end
        for (k, a) in enumerate(mg.grids)
            if a.carrier isa type
                plot(hours, a.carrier.power[:,y,s], label = string("Grids ", k))
            end
        end
        legend()
    end
end
# Statistics
function plot_metrics(metrics::Metrics; type = "hist")
    # Seaborn configuration
    Seaborn.set_theme(context="notebook", style="ticks", palette="muted", font="serif", font_scale = 1.5)

    if type == "hist"
        figure("Renewable share")
        hist(reshape(metrics.renewable_share[2:end, :], :) * 100)
        ylabel("Counts", size = "large"), yticks(size = "medium")
        xlabel("Renewable share (%)", size = "large"), xticks(size = "medium")
        figure("Cost")
        hist(reshape(metrics.eac.total[:, :], :) / 1000)
        ylabel("Counts", size = "large"), yticks(size = "medium")
        xlabel("Annual cost (k€/y)", size = "large"), xticks(size = "medium")
        if !isa(metrics.lpsp.elec, Nothing)
            figure("LPSP elec")
            hist(reshape(metrics.lpsp.elec[2:end, :], :) * 100)
            ylabel("Counts", size = "large"), yticks(size = "medium")
            xlabel("LPSP (%)", size = "large"), xticks(size = "medium")
        end
        if !isa(metrics.lpsp.heat, Nothing)
            figure("LPSP heat")
            hist(reshape(metrics.lpsp.heat[2:end, :], :) * 100)
            ylabel("Counts", size = "large"), yticks(size = "medium")
            xlabel("LPSP (%)", size = "large"), xticks(size = "medium")
        end
        if !isa(metrics.lpsp.EnergyCarrier, Nothing)
            figure("LPSP EnergyCarrier")
            hist(reshape(metrics.lpsp.EnergyCarrier[2:end, :], :) * 100)
            ylabel("Counts", size = "large"), yticks(size = "medium")
            xlabel("LPSP (%)", size = "large"), xticks(size = "medium")
        end
    elseif type == "violin"
        figure("Renewable share")
        violinplot(reshape(metrics.renewable_share[2:end, :], :) * 100)
        yticks(size = "medium")
        xlabel("Renewable share (%)", size = "large"), xticks(size = "medium")
        figure("Cost")
        violinplot(reshape(metrics.eac.total[:, :], :) / 1000)
        yticks(size = "medium")
        xlabel("Annual cost (k€/y)", size = "large"), xticks(size = "medium")
        if !isa(metrics.lpsp.elec, Nothing)
            figure("LPSP elec")
            violinplot(reshape(metrics.lpsp.elec[2:end, :], :) * 100)
            yticks(size = "medium")
            xlabel("LPSP (%)", size = "large"), xticks(size = "medium")
        end
        if !isa(metrics.lpsp.heat, Nothing)
            figure("LPSP heat")
            violinplot(reshape(metrics.lpsp.heat[2:end, :], :) * 100)
            yticks(size = "medium")
            xlabel("LPSP (%)", size = "large"), xticks(size = "medium")
        end
        if !isa(metrics.lpsp.EnergyCarrier, Nothing)
            figure("LPSP EnergyCarrier")
            violinplot(reshape(metrics.lpsp.EnergyCarrier[2:end, :], :) * 100)
            yticks(size = "medium")
            xlabel("LPSP (%)", size = "large"), xticks(size = "medium")
        end
    else
        println("Only 'hist' or 'violin' type accepted")
    end
end
