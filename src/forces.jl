"""Lenard Jones forces"""
function lennard_jones_force(r::SVector{3,Float64})
    r2 = sum(abs2, r)

    inv_r2 = 1.0 / r2
    r4 = inv_r2^2
    r8 = r4^2
    r14 = r8 * r4 * inv_r2

    return r * (r14 - 0.5*r8)
end


function calculate_force(positions, cutoff, L)

    N=length(positions)

    forces = fill(zero(SVector{3,Float64}), N)
    cutoff2=cutoff^2 # Again, use square do compare distances

    @inbounds for i in 1:N-1

        p_i=positions[i]
        for j=i+1:N

            r = minimum_image(p_i-positions[j], L)

            r2=sum(abs2,r)
            if r2 < cutoff2 # compare square distances to avoid norms
                fij=lennard_jones_force(r)

                forces[i]+=fij
                forces[j]-=fij
            end
        end
    end

    return forces
end