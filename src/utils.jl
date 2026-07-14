"""
Generate FCC lattice positions
"""
function generate_fcc_lattice(L::Float64, Ncell::Int)

    N = 4*Ncell^3
    positions = Vector{SVector{3,Float64}}(undef,N)

    idx = 1
    shift = 0.5

    for x in 1:Ncell
        for y in 1:Ncell
            for z in 1:Ncell
                positions[idx] = SVector(x-shift, y-shift, z-shift) / Ncell * L
                idx += 1

                positions[idx] =
                    SVector(x-shift, y, z) / Ncell * L
                idx += 1

                positions[idx] =
                    SVector(x, y-shift, z) / Ncell * L
                idx += 1

                positions[idx] =
                    SVector(x, y, z-shift) / Ncell * L
                idx += 1
            end
        end
    end

    println("generated ",N," particle positions")

    return positions
end




function collect_pair_distances( positions, lut, L, cutoff)

    distances = Float64[]
    cutoff2 = cutoff^2 # use square distances


    @inbounds for i in eachindex(positions)
        p_i = positions[i]
        for j in lut[i]
            r = minimum_image(p_i - positions[j], L)

            r2 = r ⋅ r
            if r2 < cutoff2
                d = sqrt(r2)
                # one entry for each particle
                push!(distances,d)
                push!(distances,d)
            end

        end
    end

    return distances
end



function getVelocities(
    T0::Float64,
    N::Int;
    seed=123
)

    Random.seed!(seed)

    velocities = Vector{SVector{3,Float64}}(undef,N)

    if T0 == 0
        fill!(velocities,zero(SVector{3,Float64}))
        return velocities
    end

    sigma = sqrt(T0)/sqrt(48)

    d = Normal(0,sigma)

    vcm = zeros(3)

    for i in 1:N
        v = SVector(rand(d), rand(d), rand(d))
        velocities[i] = v
        vcm .+= v
    end
 
    vcm ./= N
    cm = SVector{3,Float64}(vcm)

    for i in 1:N
        velocities[i]-=cm
    end

    return velocities
end



function thermostat!(Tgoal, velocities)

    T = 16 * mean(sum(abs2,v) for v in velocities)

    scale=sqrt(Tgoal/T)

    @inbounds for i in eachindex(velocities)
        velocities[i]*=scale
    end

end

