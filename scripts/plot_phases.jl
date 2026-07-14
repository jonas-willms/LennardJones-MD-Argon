using ArgonMD
using StaticArrays  
using GLMakie      
using StatsBase
###### Example script to generate a plot showing different states of argon and the neighbor distributions ######

# Investigate three different combinations of rho and temperature
rhos = [1.2, 0.8, 0.5] 
T_0s = [0.3, 0.8, 1.8]

N = 864 # Number of particles

dt = 0.032 #Time step
m = 1.0 # Mass
 
t_end = 25.0 # Total time

# Limits for potential and LUT
cutoff_Lut = 3.3
cutoff_potential = 2.5
pair_distance_limit = 2.5
n_lut_updates = 16

evaluation_steps = 30 # evaluate last evaluation_steps steps at the end of the simulation

hist_range = 0:0.05:pair_distance_limit # binning for the histogram

n_steps = Int(ceil(t_end/dt))
println("Total integration steps: ", n_steps)


if n_steps-evaluation_steps < 0 
    @error "cannot evaluate more steps than the simualtion has... $(n_steps) < $(evaluation_steps)"
end


equilibrated_positions = Vector{Vector{SVector{3,Float64}}}()

# Store histogramms in here
hs = Vector{Vector{Float64}}()

# Evolve the three param combinations
for state in 1:3
    # Pick params
    rho = rhos[state]
    T0 = T_0s[state]

    L = (N/rho)^(1/3) # Determine domain size

    println("--------------------------------")
    println("State $state: ρ=$rho T=$T0 L=$(round(L,digits=3))")

    # initial positions
    pos = generate_fcc_lattice(L, 6)


    # velocities
    vel = getVelocities(T0, length(pos))

    # initial forces
    forces = calculate_force(pos, cutoff_potential, L)

    # initial neighbour list
    lut = update_lut(pos, cutoff_Lut, L)

    state_distances = Float64[]

    # integration
    for step in 1:n_steps
        if step % n_lut_updates == 0
            lut = update_lut(pos, cutoff_Lut, L)
        end

        # Equilibration thermostat only for
        # the first half of the simulation
        # every 10 steps
        # and T0 > 0
        if step < n_steps÷2 && step % 10 == 0 && T0 > 0
            thermostat!(T0, vel)
        end

        # Evolve the system one dt
        step_lut!(pos, vel, forces, lut, dt, L, cutoff_potential, m)

        # collect final statistics
        if step > n_steps-evaluation_steps
            append!(state_distances, collect_pair_distances(pos, lut, L, pair_distance_limit))
        end

    end


    println("finished state $state, collected ",length(state_distances)," distances")

    # save positions
    push!(equilibrated_positions, copy(pos))

    # histogram
    hist = fit(Histogram, state_distances, hist_range)
    h = hist # renember this as baseline for the plotting edges

    push!(hs, hist.weights ./ (N *  evaluation_steps))

end


println("Simulation complete")

# ============================================================
# Plotting
# ============================================================

GLMakie.activate!(inline=true)

set_theme!(
    Theme(
        fontsize = 18,
        Axis = (
            xlabelsize = 22,
            ylabelsize = 22,
            titlesize = 24,
            xticklabelsize = 16,
            yticklabelsize = 16,
        )
    )
)


fig = Figure(
    size = (1000, 1000),
    fontsize = 18
)


phase_names = [
    "solid",
    "liquid",
    "gas"
]


phase_colors = [
    :red,
    :green,
    :blue
]



# ------------------------------------------------------------
# Pair distribution histograms
# ------------------------------------------------------------
xs = [x + (hist_range[end] - hist_range[1])/length(hist_range) for x in hist_range[1:end-1]]

for i in 1:3

    ax = Axis(fig[i,1:2],
        xlabel = L"particle pair distance $r$",
        ylabel = L"#Neighbors $n(r)$",
        title = L"Density $\rho=%$(rhos[i])$,\quad Temperature $T=%$(T_0s[i])$"
    )


    barplot!(ax,
        xs,
        hs[i],
        color = phase_colors[i],
        strokecolor = :black,
        strokewidth = 1,
        alpha = 0.85
    )


    text!(
        ax,
        0.25,
        maximum(hs[i])*0.7,
        text = phase_names[i],
        fontsize = 24,
        color = :black
    )

end



# ------------------------------------------------------------
# Particle snapshots
# ------------------------------------------------------------
a = 5.25 - 2.5
b = 5.25 + 2.5


for i in 1:3
    ax = Axis(fig[i,3],
        xlabel = L"$x$",
        ylabel = L"$y$",
        aspect = DataAspect(),
        title = L"\rho=%$(rhos[i]),\quad T=%$(T_0s[i])",
        limits = (a,b,a,b),
    )


    pts = [p for p in equilibrated_positions[i] if ( a < p[1] < b && a < p[2] < b && a < p[3] < b)]

    scatter!(
        ax,
        first.(pts),
        getindex.(pts,2),
        markersize = 8,
        color = [
            (
                phase_colors[i],
                clamp((p[3]-2)/5,0,1)
            )
            for p in pts
        ]
    )

end


#save the image
save("states.png", fig, update = true, px_per_unit = 2)