using JuMP
import Ipopt


# Load data and model
include("data.jl")
the_model = Model(Ipopt.Optimizer)


# Variables
@variable(the_model, 0 <= Gp[i = 1:n_gen] <= G_cap[i])
@variable(the_model, -G_cap[i]*0.003 <= Gq[i = 1:n_gen] <= G_cap[i]*0.003)
@variable(the_model, C[i = 1:n_cust] >= C_dem[i])
@variable(the_model, -pi <= theta[i = 1:n_nodes] <= pi)
@variable(the_model, 0.98 <= v[i = 1:n_nodes] <= 1.02)
# @variable(the_model, p[i = 1:n_flows])
# @variable(the_model, q[i = 1:n_flows])
@variable(the_model, p[i = 1:n_nodes, j = 1:n_nodes])
@variable(the_model, q[i = 1:n_nodes, j = 1:n_nodes])


# Objective function
@NLobjective(the_model, Min, sum((Gp[i]*G_cost[i]) for i in 1:n_gen))


# Netflow constraints from node k -> l
for k in 1:n_nodes
    for l in 1:n_nodes
        # Mirroring the initial b and g arrays
        if l > k
            b[l,k] = b[k,l]
            g[l,k] = g[k,l]
        end
        @NLconstraint(the_model, p[k,l] ==  (v[k]^2)*g[k,l] - v[k]*v[l]*g[k,l]*cos(theta[k]-theta[l]) - v[k]*v[l]*b[k,l]*sin(theta[k]-theta[l]))
        @NLconstraint(the_model, q[k,l] == -(v[k]^2)*b[k,l] + v[k]*v[l]*b[k,l]*cos(theta[k]-theta[l]) - v[k]*v[l]*g[k,l]*sin(theta[k]-theta[l]))
    end
end


# Flow balance constraints for every node: generated power + incoming power = absorbed power
# Active power
@NLconstraint(the_model, np1, 0 == C[1] + p[1,2] + p[1,11])
@NLconstraint(the_model, np2, Gp[1] + Gp[2] + Gp[3] == p[2,1] + p[2,11] + p[2,3])
@NLconstraint(the_model, np3, Gp[4] == p[3,2] + p[3,4] + p[3,9])
@NLconstraint(the_model, np4, Gp[5] == C[2] + p[4,3] + p[4,5])
@NLconstraint(the_model, np5, Gp[6] ==  p[5,4] + p[5,6] + p[5,8])
@NLconstraint(the_model, np6, 0 == C[3] + p[6,5] + p[6,7])
@NLconstraint(the_model, np7, Gp[7] ==  p[7,6] + p[7,8] + p[7,9])
@NLconstraint(the_model, np8, 0 == C[4] + p[8,5] + p[8,7] + p[8,9])
@NLconstraint(the_model, np9, Gp[8] + Gp[9] == C[5] + p[9,3] + p[9,7] + p[9,8] + p[9,10])
@NLconstraint(the_model, np10, 0 == C[6] + p[10,9] + p[10,11])
@NLconstraint(the_model, np11, 0  == C[7] + p[11,1] + p[11,2] + p[11,10])

# Reactive power
@NLconstraint(the_model, nq1, 0 == q[1,2] + q[1,11])
@NLconstraint(the_model, nq2, Gq[1] + Gq[2] + Gq[3] == q[2,1] + q[2,11] + q[2,3])
@NLconstraint(the_model, nq3, Gq[4] == q[3,2] + q[3,4] + q[3,9])
@NLconstraint(the_model, nq4, Gq[5] == q[4,3] + q[4,5])
@NLconstraint(the_model, nq5, Gq[6] == q[5,4] + q[5,6] + q[5,8])
@NLconstraint(the_model, nq6, 0 == q[6,5] + q[6,7])
@NLconstraint(the_model, nq7, Gq[7] == q[7,6] + q[7,8] + q[7,9])
@NLconstraint(the_model, nq8, 0 == q[8,5] + q[8,7] + q[8,9])
@NLconstraint(the_model, nq9, Gq[8] + Gq[9] == q[9,3] + q[9,7] + q[9,8] + q[9,10])
@NLconstraint(the_model, nq10, 0 == q[10,9] + q[10,11])
@NLconstraint(the_model, nq11, 0  == q[11,1] + q[11,2] + q[11,10])


# Print and optimize model
println(the_model)
optimize!(the_model)


# Printing some of the results for further analysis
println("") # Printing white line after solver output, before printing
println("Termination statue: ", JuMP.termination_status(the_model))
println("Optimal objective function value: ", JuMP.objective_value(the_model))
println("Generators active power: ", JuMP.value.(Gp[i] for i in 1:n_gen))
println("Generators reactive power: ", JuMP.value.(Gq[i] for i in 1:n_gen))
println("Customers active power: ", JuMP.value.(C[i] for i in 1:n_cust))

println("Pkl: ", JuMP.value.(p[1,i] for i in 1:n_nodes))
println("Pkl: ", JuMP.value.(p[2,i] for i in 1:n_nodes))
println("Pkl: ", JuMP.value.(p[3,i] for i in 1:n_nodes))
println("Pkl: ", JuMP.value.(p[4,i] for i in 1:n_nodes))
println("Pkl: ", JuMP.value.(p[5,i] for i in 1:n_nodes))
println("Pkl: ", JuMP.value.(p[6,i] for i in 1:n_nodes))
println("Pkl: ", JuMP.value.(p[7,i] for i in 1:n_nodes))
println("Pkl: ", JuMP.value.(p[8,i] for i in 1:n_nodes))
println("Pkl: ", JuMP.value.(p[9,i] for i in 1:n_nodes))
println("Pkl: ", JuMP.value.(p[10,i] for i in 1:n_nodes))
println("Pkl: ", JuMP.value.(p[11,i] for i in 1:n_nodes))

println("qkl: ", JuMP.value.(q[1,i] for i in 1:n_nodes))
println("qkl: ", JuMP.value.(q[2,i] for i in 1:n_nodes))
println("qkl: ", JuMP.value.(q[3,i] for i in 1:n_nodes))
println("qkl: ", JuMP.value.(q[4,i] for i in 1:n_nodes))
println("qkl: ", JuMP.value.(q[5,i] for i in 1:n_nodes))
println("qkl: ", JuMP.value.(q[6,i] for i in 1:n_nodes))
println("qkl: ", JuMP.value.(q[7,i] for i in 1:n_nodes))
println("qkl: ", JuMP.value.(q[8,i] for i in 1:n_nodes))
println("qkl: ", JuMP.value.(q[9,i] for i in 1:n_nodes))
println("qkl: ", JuMP.value.(q[10,i] for i in 1:n_nodes))
println("qkl: ", JuMP.value.(q[11,i] for i in 1:n_nodes))

println("Dual variables/Lagrange multipliers corresponding to some of the constraints: ")
# println(JuMP.dual.(SOS_constr))
println(JuMP.dual.(JuMP.UpperBoundRef.(Gp)))
