# Load packages
using Genesys, JLD, Dates

#=
First, you have to define your own designer and controller type. Make sure that the genesys "supertype" is added.
=#

# Define your own designer...
mutable struct bar <: Genesys.AbstractDesigner
    decisions::NamedTuple
    bar() = new()
end

# ...and controller
mutable struct foo <: Genesys.AbstractController
    decisions::NamedTuple
    foo() = new()
end

# Define the "offline" functions for both the designer and controller
function Genesys.initialize_controller!(mg::Microgrid, controller::foo, ω::Scenarios)
    # Preallocate
    Genesys.preallocate!(mg, controller)
    return controller
end

function Genesys.initialize_designer!(mg::Microgrid, designer::bar, ω::Scenarios)
    # Preallocate
    Genesys.preallocate!(mg, designer)
    return designer
end

# Define the "online" functions for both the designer and controller
function Genesys.compute_operation_decisions!(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::foo)
    return controller
end

function Genesys.compute_investment_decisions!(y::Int64, s::Int64, mg::Microgrid, designer::bar)
    return designer
end

#=
Let's simulate the microgrid with the dummies controller and designer...
=#

# Parameters of the simulation
const nh, ny, ns = 8760, 2, 1

# Load input data
data = load(joinpath("data","ausgrid_5_twostage.jld"))

# Initialize the microgrid
microgrid = Microgrid(parameters = GlobalParameters(nh, ny, ns))

# Add the equipment to the microgrid
add!(microgrid, Demand(carrier = Electricity()), Solar(), Liion(), Grid(carrier = Electricity()))

# Initialize scenarios
ω_d, ω_a = Scenarios(microgrid, data["ω_optim"]), Scenarios(microgrid, data["ω_simu"])

# Initialize designer
designer = initialize_designer!(microgrid, bar(), ω_d)

# Initialize controller
controller = initialize_controller!(microgrid, foo(), ω_d)

# Assessment
simulate!(microgrid, controller, designer, ω_a)
