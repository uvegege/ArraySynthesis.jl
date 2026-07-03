
# U: [3 x N] direction matrix with columns (u, v, w)
function steering_matrix(a::ArrayGeometry, U, k = 2 * π)
    ψ = U' * a.positions
    return cis.(k*ψ)
end

function steering_matrix(a::SymmetricArray, U, k = 2 * π)
    ψ = k .* (U' * a.positions)
    A_cos = 2 .* cos.(ψ)
    A_sin = -2 .* sin.(ψ)
    return A_cos, A_sin
end
