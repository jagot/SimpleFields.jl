module SimpleFields
using FixedSizeArrays
import Base: call

abstract Field{U<:AbstractFloat}

# λ_SI in meters
function fundamental{U<:AbstractFloat}(λ_SI::U)
    λ = λ_SI/5.2917721067e-11
    T = λ/137.035999
    ω = 2π/T
    λ,T,ω
end

# I_SI in W/cm²
intensity{U<:AbstractFloat}(I_SI::U) = I_SI/3.5094452e16

# Gaussian envelope for the amplitude, fwhm of the intensity envelope
# in cycles.
gaussian{U<:AbstractFloat}(t::U, fwhm::U) = exp(-t.^2/(fwhm/(2*√(log(2)))))

type LinearField{U<:AbstractFloat} <: Field{U}
    λ::U
    T::U
    ω::U
    tmax::U
    Ez::Function
end

call{U<:AbstractFloat}(field::LinearField, t::U) = field.Ez(t)

type TransverseField{U<:AbstractFloat} <: Field{U}
    λ::U
    T::U
    ω::U
    tmax::U
    Ez::Function
    Ex::Function
end

call{U<:AbstractFloat}(field::TransverseField, t::U) = Vec{2,Float64}(field.Ez(t),field.Ex(t))

call{U<:AbstractFloat}(field::Field{U}, t::AbstractVector{U}) = map(field, t)

# λ_SI in meters, I_SI in W/cm², tmax in cycles, fwhm in seconds, cep
# in units of π. The pulse will be centered at tmax/2.
function pulse{U<:AbstractFloat}(λ_SI::U, I_SI::U,
                                 tmax::U, fwhm::U,
                                 env = gaussian,
                                 cep = 0.0)
    λ,T,ω = fundamental(λ_SI)
    I = intensity(I_SI)
    E = sqrt(I)
    fwhm /= 2.41888430e-17T
    LinearField(λ, T, ω, tmax, t -> E*gaussian(t-tmax/2, fwhm)*sin(2π*(t-tmax/2) + cep*π))
end

export Field, pulse, call
end # module
