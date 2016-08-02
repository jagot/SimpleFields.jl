module SimpleFields
using FixedSizeArrays
import Base: eltype, call, zero, (+)

abstract Field{U<:AbstractFloat}
eltype{U<:AbstractFloat}(::Type{Field{U}}) = U

type CompositeField{U<:AbstractFloat} <: Field{U}
    components::Vector{Field}
    λ::U
    T::U
    ω::U
    tmax::U
    zero
end

function call{U<:AbstractFloat}(field::CompositeField{U}, t::U)
    f = field.zero
    for c in field.components
        f += c(t)
    end
    f
end

function (+){F1<:Field,F2<:Field}(a::F1, b::F2)
    zero(F1) == zero(F2) || error("Can only add fields of same type")
    a.λ == b.λ || error("Can only add fields of the same fundamental wavelength")
    CompositeField{eltype(a)}([a, b], a.λ, a.T, a.ω, max(a.tmax, b.tmax), zero(F1))
end

function (+){F<:Field}(a::CompositeField, b::F)
    a.zero == zero(typeof(b)) || error("Can only add fields of same type")
    a.λ == b.λ || error("Can only add fields of the same fundamental wavelength")
    CompositeField{eltype(a)}([a.components..., b], a.λ, a.T, a.ω, max(a.tmax, b.tmax), a.zero)
end

type DelayedField{U<:AbstractFloat,F<:Field} <: Field{U}
    f::F
    tau::U
    λ::U
    T::U
    ω::U
    tmax::U
end
eltype{U<:AbstractFloat,F<:Field}(::Type{DelayedField{U,F}}) = eltype(F)
call{U<:AbstractFloat}(field::DelayedField{U}, t::U) = field.f(t-field.tau)
zero{U<:AbstractFloat,F<:Field}(::Type{DelayedField{U,F}}) = zero(F)

delay{U<:AbstractFloat, F<:Field}(f::F, tau::U) = DelayedField{U,F}(f, tau, f.λ, f.T, f.ω, f.tmax+tau)

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
gaussian{U<:AbstractFloat}(t::U, fwhm::U) = exp(-t.^2/(2(fwhm/(2*√(log(2))))^2))

box{U<:AbstractFloat}(t::U, tmax::U, c::Bool) = !c || (t >= 0 && t <= tmax) ? one(U) : zero(U)

type LinearField{U<:AbstractFloat} <: Field{U}
    λ::U
    T::U
    ω::U
    tmax::U
    Ez::Function
    vanish::Bool
end

call{U<:AbstractFloat}(field::LinearField{U}, t::U) = field.Ez(t)*box(t,field.tmax,field.vanish)
zero{U<:AbstractFloat}(::Type{LinearField{U}}) = zero(U)

type TransverseField{U<:AbstractFloat} <: Field{U}
    λ::U
    T::U
    ω::U
    tmax::U
    Ez::Function
    Ex::Function
    vanish::Bool
end

call{U<:AbstractFloat}(field::TransverseField{U}, t::U) = Vec{2,U}(field.Ez(t),field.Ex(t))*box(t,field.tmax,field.vanish)
zero{U<:AbstractFloat}(::Type{TransverseField{U}}) = Vec{2,U}(zero(U), zero(U))

call{U<:AbstractFloat}(field::Field{U}, t::AbstractVector{U}) = map(field, t)

function gdd_params{U<:AbstractFloat}(λ_SI::U, τ₀::U, η::U = 0.0;
                                      gdd_phase = false)
    γ = τ₀^2/(8log(2))
    γ² = γ^2
    η² = η^2

    A = √(γ/√(γ² + η²))
    gdd_phase && (A *= exp(im*atan2(-η,γ)/2))

    a = 1./4*(γ/(γ² + η²))
    b = 1./2*(η/(γ² + η²))
    A,a,b
end

# λ_SI in meters, I_SI in W/cm², tmax in cycles, (Fourier-limited)
# fwhm in seconds, cep in units of π, gdd in seconds squared. The
# pulse will be centered at tmax/2. If vanish is true, the pulse will
# identically vanish outside the interval [0,tmax]
function pulse{U<:AbstractFloat}(λ_SI::U, I_SI::U,
                                 tmax::U, fwhm::U,
                                 q::U = 1.0,
                                 cep::U = 0.0, gdd::U = 0.0;
                                 gdd_phase = false,
                                 vanish = true)
    λ,T,ω = fundamental(λ_SI)
    T2 = (2.41888430e-17T)^2
    I = intensity(I_SI)
    E₀ = √(I)

    A,a,b = gdd_params(λ_SI, fwhm, gdd)

    A *= E₀

    # Transfer to cycles as time base
    a *= T2
    b *= T2

    φ₀ = cep*π

    E = t -> A*exp(im*2π*q*(t-tmax/2) + im*q*φ₀)*exp(-(a-im*b)*(t-tmax/2)^2)

    LinearField(λ, T, ω, tmax, t -> real(E(t)), vanish)
end

export Field, CompositeField, fundamental, delay, gdd_params, pulse, eltype, call, (+)
end # module
