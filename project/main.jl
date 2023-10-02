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
@variable(the_model, p[i = 1:n_flows])
@variable(the_model, q[i = 1:n_flows])


# Objective function
@NLobjective(the_model, Min, sum((Gp[i]*G_cost[i]) for i in 1:n_gen))


# Netflow constraints from node k -> l
# Active power
@NLconstraint(the_model, p1, p[1] == g_kl[1]*(v[2]^2-v[1]^2) - 2*b_kl[1]*v[2]*v[1]*sin(theta[2]-theta[1])) # 2 -> 1
@NLconstraint(the_model, p2, p[2] == g_kl[2]*(v[11]^2-v[1]^2) - 2*b_kl[2]*v[11]*v[1]*sin(theta[11]-theta[1])) # 11 -> 1
@NLconstraint(the_model, p3, p[3] == g_kl[3]*(v[3]^2-v[2]^2) - 2*b_kl[3]*v[3]*v[2]*sin(theta[3]-theta[2])) # 3 -> 2
@NLconstraint(the_model, p4, p[4] == g_kl[4]*(v[11]^2-v[2]^2) - 2*b_kl[4]*v[11]*v[2]*sin(theta[11]-theta[2])) # 11 -> 2
@NLconstraint(the_model, p5, p[5] == g_kl[5]*(v[4]^2-v[3]^2) - 2*b_kl[5]*v[4]*v[3]*sin(theta[4]-theta[3])) # 4 -> 3
@NLconstraint(the_model, p6, p[6] == g_kl[6]*(v[9]^2-v[3]^2) - 2*b_kl[6]*v[9]*v[3]*sin(theta[9]-theta[3])) # 9 -> 3
@NLconstraint(the_model, p7, p[7] == g_kl[7]*(v[5]^2-v[4]^2) - 2*b_kl[7]*v[5]*v[4]*sin(theta[5]-theta[4])) # 5 -> 4
@NLconstraint(the_model, p8, p[8] == g_kl[8]*(v[6]^2-v[5]^2) - 2*b_kl[8]*v[6]*v[5]*sin(theta[6]-theta[5])) # 6 -> 5
@NLconstraint(the_model, p9, p[9] == g_kl[9]*(v[8]^2-v[5]^2) - 2*b_kl[9]*v[8]*v[5]*sin(theta[8]-theta[5])) # 8 -> 5
@NLconstraint(the_model, p10, p[10] == g_kl[10]*(v[7]^2-v[6]^2) - 2*b_kl[10]*v[7]*v[6]*sin(theta[7]-theta[6])) # 7 -> 6
@NLconstraint(the_model, p11, p[11] == g_kl[11]*(v[8]^2-v[7]^2) - 2*b_kl[11]*v[8]*v[7]*sin(theta[8]-theta[7])) # 8 -> 7 
@NLconstraint(the_model, p12, p[12] == g_kl[12]*(v[9]^2-v[7]^2) - 2*b_kl[12]*v[9]*v[7]*sin(theta[9]-theta[7])) # 9 -> 7
@NLconstraint(the_model, p13, p[13] == g_kl[13]*(v[9]^2-v[8]^2) - 2*b_kl[13]*v[9]*v[8]*sin(theta[9]-theta[8])) # 9 -> 8
@NLconstraint(the_model, p14, p[14] == g_kl[14]*(v[10]^2-v[9]^2) - 2*b_kl[14]*v[10]*v[9]*sin(theta[10]-theta[9])) # 10 -> 9
@NLconstraint(the_model, p15, p[15] == g_kl[15]*(v[11]^2-v[10]^2) - 2*b_kl[15]*v[11]*v[10]*sin(theta[11]-theta[10])) # 11 -> 10

# Reactive power
@NLconstraint(the_model, q1, q[1] == b_kl[1]*(v[1]^2-v[2]^2) - 2*g_kl[1]*v[2]*v[1]*sin(theta[2]-theta[1])) # 2 -> 1
@NLconstraint(the_model, q2, q[2] == b_kl[2]*(v[1]^2-v[11]^2) - 2*g_kl[2]*v[11]*v[1]*sin(theta[11]-theta[1])) # 11 -> 1
@NLconstraint(the_model, q3, q[3] == b_kl[3]*(v[2]^2-v[3]^2) - 2*g_kl[3]*v[3]*v[2]*sin(theta[3]-theta[2])) # 3 -> 2
@NLconstraint(the_model, q4, q[4] == b_kl[4]*(v[2]^2-v[11]^2) - 2*g_kl[4]*v[11]*v[2]*sin(theta[11]-theta[2])) # 11 -> 2
@NLconstraint(the_model, q5, q[5] == b_kl[5]*(v[3]^2-v[4]^2) - 2*g_kl[5]*v[4]*v[3]*sin(theta[4]-theta[3])) # 4 -> 3
@NLconstraint(the_model, q6, q[6] == b_kl[6]*(v[3]^2-v[9]^2) - 2*g_kl[6]*v[9]*v[3]*sin(theta[9]-theta[3])) # 9 -> 3
@NLconstraint(the_model, q7, q[7] == b_kl[7]*(v[4]^2-v[5]^2) - 2*g_kl[7]*v[5]*v[4]*sin(theta[5]-theta[4])) # 5 -> 4
@NLconstraint(the_model, q8, q[8] == b_kl[8]*(v[5]^2-v[6]^2) - 2*g_kl[8]*v[6]*v[5]*sin(theta[6]-theta[5])) # 6 -> 5
@NLconstraint(the_model, q9, q[9] == b_kl[9]*(v[5]^2-v[8]^2) - 2*g_kl[9]*v[8]*v[5]*sin(theta[8]-theta[5])) # 8 -> 5
@NLconstraint(the_model, q10, q[10] == b_kl[10]*(v[6]^2-v[7]^2) - 2*g_kl[10]*v[7]*v[6]*sin(theta[7]-theta[6])) # 7 -> 6
@NLconstraint(the_model, q11, q[11] == b_kl[11]*(v[7]^2-v[8]^2) - 2*g_kl[11]*v[8]*v[7]*sin(theta[8]-theta[7])) # 8 -> 7 
@NLconstraint(the_model, q12, q[12] == b_kl[12]*(v[7]^2-v[9]^2) - 2*g_kl[12]*v[9]*v[7]*sin(theta[9]-theta[7])) # 9 -> 7
@NLconstraint(the_model, q13, q[13] == b_kl[13]*(v[8]^2-v[9]^2) - 2*g_kl[13]*v[9]*v[8]*sin(theta[9]-theta[8])) # 9 -> 8
@NLconstraint(the_model, q14, q[14] == b_kl[14]*(v[9]^2-v[10]^2) - 2*g_kl[14]*v[10]*v[9]*sin(theta[10]-theta[9])) # 10 -> 9
@NLconstraint(the_model, q15, q[15] == b_kl[15]*(v[10]^2-v[11]^2) - 2*g_kl[15]*v[11]*v[10]*sin(theta[11]-theta[10])) # 11 -> 10


# Flow balance constraints of every node: generated power - absorbed power = netinflow - netoutflow
# Active power
@NLconstraint(the_model, np1, C[1] == p[1] + p[2])
@NLconstraint(the_model, np2, Gp[1] + Gp[2] + Gp[3] == p[1] -p[3] -p[4])
@NLconstraint(the_model, np3, Gp[4] == p[3] -p[5] -p[6])
@NLconstraint(the_model, np4, Gp[5] - C[2] == p[5] -p[7])
@NLconstraint(the_model, np5, Gp[6] == p[7] -p[8] -p[9])
@NLconstraint(the_model, np6, C[3] == -p[8] + p[10])
@NLconstraint(the_model, np7, Gp[7] == p[10] -p[11] -p[12])
@NLconstraint(the_model, np8, C[4] == -p[9] -p[11] + p[13])
@NLconstraint(the_model, np9, Gp[8] + Gp[9] - C[5] == p[6] + p[12] + p[13] -p[14])
@NLconstraint(the_model, np10, C[6] == -p[14] + p[15])
@NLconstraint(the_model, np11, C[7] == -p[2] -p[4] -p[15])

# Reactive power
@NLconstraint(the_model, nq1, 0 == q[1] + q[2])
@NLconstraint(the_model, nq2, Gq[1] + Gq[2] + Gq[3] == q[1] -q[3] -q[4])
@NLconstraint(the_model, nq3, Gq[4] == q[3] -q[5] -q[6])
@NLconstraint(the_model, nq4, Gq[5] == q[5] -q[7])
@NLconstraint(the_model, nq5, Gq[6] == q[7] -q[8] -q[9])
@NLconstraint(the_model, nq6, 0 == -q[8] + q[10])
@NLconstraint(the_model, nq7, Gq[7] == q[10] -q[11] -q[12])
@NLconstraint(the_model, nq8, 0 == -q[9] -q[11] + q[13])
@NLconstraint(the_model, nq9, Gq[8] + Gq[9] == q[6] + q[12] + q[13] -q[14])
@NLconstraint(the_model, nq10, 0 == -q[14] + q[15])
@NLconstraint(the_model, nq11, 0 == -q[2] -q[4] -q[15])


# Print and optimize model
println(the_model)
optimize!(the_model)


# Printing some of the results for further analysis
println("") # Printing white line after solver output, before printing
println("Termination statue: ", JuMP.termination_status(the_model))
println("Optimal objective function value: ", JuMP.objective_value(the_model))
println("Generators activ power: ", JuMP.value.(Gp[i] for i in 1:n_gen))
println("Customers active power: ", JuMP.value.(C[i] for i in 1:n_cust))

println("Dual variables/Lagrange multipliers corresponding to some of the constraints: ")
# println(JuMP.dual.(SOS_constr))
println(JuMP.dual.(JuMP.UpperBoundRef.(Gp)))
