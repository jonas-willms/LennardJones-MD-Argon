
function update_lut(positions, cutoff, L)

    N=length(positions)

    # Re-allocating the memory as the size of the entries may vary
    lut = [Int[] for _ in 1:N]

    cutoff2=cutoff^2 # Squared distances for comaprisons

    @inbounds for i in 1:N-1
        p_i=positions[i]

        for j=i+1:N
            r = minimum_image(p_i-positions[j], L)

            if sum(abs2,r)<cutoff2
                push!(lut[i],j)
            end
        end
    end

    return lut
end

function minimum_image(r::SVector{3,Float64}, L::Float64)
    return r - L * SVector(
        round(r[1]/L),
        round(r[2]/L),
        round(r[3]/L)
    )
end



# ============================================================
# Integration step
# ============================================================
"""Use a velocity Verlet step with a lookup taple to evolve the system"""
function step_lut!(positions, velocities, forces, lut, dt, L, cutoff, m)

    N=length(positions)

    # velocity half step + move
    @inbounds for i in 1:N
        velocities[i]+=forces[i]*dt/(2m)
        positions[i]+=velocities[i]*dt
        forces[i]=zero(SVector{3,Float64})
    end

    cutoff2=cutoff^2

    # forces
    @inbounds for i in 1:N
        p_i=positions[i]
        for j in lut[i]
            r = minimum_image(p_i-positions[j], L)

            r2=sum(abs2,r)

            if r2 < cutoff2
                fij=lennard_jones_force(r)
                forces[i]+=fij
                forces[j]-=fij
            end
        end
    end

    # second half step + PBC
    @inbounds for i in 1:N
        velocities[i]+=forces[i]*dt/(2m)
        p=positions[i]

        positions[i]=SVector( mod(p[1],L), mod(p[2],L), mod(p[3],L)) # enforce periodic boundary conditions
    end

    return nothing
end