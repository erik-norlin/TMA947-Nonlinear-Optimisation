using JuMP
import Ipopt

# Import data from the data file
include("intro_data.jl")

# Create the model object
the_model = Model(Ipopt.Optimizer)

# Create (one set of) variables, and their lower and upper bounds
@variable(the_model, lb[i] <= x[i = 1:n_vars] <= ub[i])

# Create the nonlinear objective \sum_i (x_i - 1)^2, which we want to minimize
@NLobjective(the_model, Min, sum((x[i] - 1)^2 for i in 1:n_vars))

# Create the nonlinear constraint that the sum of all variables squared should
# be less than or equal to a given constant. The constrain is given the name
# "SOS_constr", which can later be referenced in the code
@NLconstraint(the_model, SOS_constr, sum(x[i]^2 for i in 1:n_vars) <= sum_bound)

# Print the optimzation problem in the terminal
println(the_model)

# Solve the optimization problem
optimize!(the_model)

# Printing some of the results for further analysis
println("") # Printing white line after solver output, before printing
println("Termination statue: ", JuMP.termination_status(the_model))
println("Optimal(?) objective function value: ", JuMP.objective_value(the_model))
println("Optimal(?) point: ", JuMP.value.(x))
println("Dual variables/Lagrange multipliers corresponding to some of the constraints: ")
println(JuMP.dual.(SOS_constr))
println(JuMP.dual.(JuMP.UpperBoundRef.(x)))
