module SimpleFields
using FixedSizeArrays
import Base: eltype, call, zero, (+)

abstract Field{U<:AbstractFloat}
eltype{U<:AbstractFloat}(f::Field{U}) = U

type CompositeField{U<:AbstractFloat,F<:Field} <: Field{U}
    components::Vector{F}
end

function call{U<:AbstractFloat,F<:Field}(field::CompositeField{U,F}, t::U)
    f = zero(F)
    for c in field.components
        f += c(t)
    end
    f
end

(+){F<:Field}(a::F, b::F) = CompositeField{eltype(a),F}([a, b])

function (+){F<:Field}(a::CompositeField, b::F)
    eltype(a.components) == typeof(b) || error("Can only add fields of same type")
    CompositeField{eltype(a),F}([a.components..., b])
end

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

call{U<:AbstractFloat}(field::LinearField{U}, t::U) = field.Ez(t)
zero{U<:AbstractFloat}(::Type{LinearField{U}}) = zero(U)

type TransverseField{U<:AbstractFloat} <: Field{U}
    λ::U
    T::U
    ω::U
    tmax::U
    Ez::Function
    Ex::Function
end

call{U<:AbstractFloat}(field::TransverseField{U}, t::U) = Vec{2,U}(field.Ez(t),field.Ex(t))
zero{U<:AbstractFloat}(::Type{TransverseField{U}}) = Vec{2,U}(zero(U), zero(U))

call{U<:AbstractFloat}(field::Field{U}, t::AbstractVector{U}) = map(field, t)

# λ_SI in meters, I_SI in W/cm², tmax in cycles, fwhm in seconds, cep
# in units of π. The pulse will be centered at tmax/2.
function pulse{U<:AbstractFloat}(λ_SI::U, I_SI::U,
                                 tmax::U, fwhm::U,
                                 q::U = 1.0,
                                 env = gaussian,
                                 cep = 0.0)
    λ,T,ω = fundamental(λ_SI)
    I = intensity(I_SI)
    E = sqrt(I)
    fwhm /= 2.41888430e-17T
    LinearField(λ, T, ω, tmax, t -> E*env(t-tmax/2, fwhm)*sin(2π*q*(t-tmax/2) + q*cep*π))
end

export Field, CompositeField, pulse, eltype, call, (+)
end # module
