function theta_ramp(θ1, y1, θ2, y2)
    θ2 == θ1 && error("θ1 and θ2 must be different.")
    return p -> begin
        θ = p.θ
        t = (θ - θ1) / (θ2 - θ1)
        return (1 - t) * y1 + t * y2
    end
end

function u_ramp(u1, y1, u2, y2)
    u2 == u1 && error("u1 and u2 must be different.")
    return p -> begin
        u = p.u
        t = (u - u1) / (u2 - u1)
        return (1 - t) * y1 + t * y2
    end
end

function v_ramp(v1, y1, v2, y2)
    v2 == v1 && error("v1 and v2 must be different.")
    return p -> begin
        t = (p.v - v1) / (v2 - v1)
        return (1 - t) * y1 + t * y2
    end
end

function csc_values(region::Region; offset = 0.1, step = 0.05)
    return [csc(offset + step * i) for i in eachindex(region.points)]
end

# Examples
#shaped_beam(region, 1.0)
#shaped_beam(region, u_ramp(0.1, 1.0, 0.5, 0.4))
#shaped_beam(region, values)
#
#sidelobes(region, -25dB)
#sidelobes(region, u_ramp(-0.3, -25dB, 0.0, -15dB))