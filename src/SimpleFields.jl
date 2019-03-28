module SimpleFields
using StaticArrays
import Base: eltype, zero, (+)

abstract type Field{U<:AbstractFloat} end
eltype(::Type{Field{U}}) where {U<:AbstractFloat} = U

struct CompositeField{U<:AbstractFloat} <: Field{U}
    components::Vector{Field}
    λ::U
    T::U
    ω::U
    tmax::U
    zero
end

function (field::CompositeField{U})(t::U) where {U<:AbstractFloat}
    f = field.zero
    for c in field.components
        f += c(t)
    end
    f
end

function (+)(a::F1, b::F2) where {F1<:Field,F2<:Field}
    zero(F1) == zero(F2) || error("Can only add fields of same type")
    a.λ == b.λ || error("Can only add fields of the same fundamental wavelength")
    CompositeField{eltype(a)}([a, b], a.λ, a.T, a.ω, max(a.tmax, b.tmax), zero(F1))
end

function (+)(a::CompositeField, b::F) where {F<:Field}
    a.zero == zero(typeof(b)) || error("Can only add fields of same type")
    a.λ == b.λ || error("Can only add fields of the same fundamental wavelength")
    CompositeField{eltype(a)}([a.components..., b], a.λ, a.T, a.ω, max(a.tmax, b.tmax), a.zero)
end

struct DelayedField{U<:AbstractFloat,F<:Field} <: Field{U}
    f::F
    tau::U
    λ::U
    T::U
    ω::U
    tmax::U
end
eltype(::Type{DelayedField{U,F}}) where {U<:AbstractFloat,F<:Field} = eltype(F)
(field::DelayedField{U})(t::U) where {U<:AbstractFloat} = field.f(t-field.tau)
zero(::Type{DelayedField{U,F}}) where {U<:AbstractFloat,F<:Field} = zero(F)

delay(f::F, tau::U) where {U<:AbstractFloat, F<:Field} = DelayedField{U,F}(f, tau, f.λ, f.T, f.ω, f.tmax+tau)

# λ_SI in meters
function fundamental(λ_SI::U) where {U<:AbstractFloat}
    λ = λ_SI/U(5.2917721067e-11)
    T = λ/U(137.035999)
    ω = U(2π/T)
    λ,T,ω
end

# I_SI in W/cm²
intensity(I_SI::U) where {U<:AbstractFloat} = I_SI/U(3.5094452e16)

# Gaussian envelope for the amplitude, fwhm of the intensity envelope
# in cycles.
gaussian(t::U, fwhm::U) where {U<:AbstractFloat} = exp(-t.^2/(2(fwhm/U(2*√(log(2))))^2))

function top_hat(t, ramp, tmax)
    ((t.>=ramp) .* (t.<tmax-ramp)) +
        ((t.>=0) .* (t.<ramp)).*(t/ramp) +
        (((t.>=tmax-ramp) .* (t.<tmax))).*((tmax-t)/ramp)
end

function top_hat_trunc(t, ramp, tmax)
    ((t.>=ramp) .* (t.<tmax)) +
        ((t.>=0) .* (t.<ramp)).*(t/ramp)
end

box(t::U, tmax::U, c::Bool) where {U<:AbstractFloat} = !c || (t >= 0 && t <= tmax) ? one(U) : zero(U)

struct LinearField{U<:AbstractFloat} <: Field{U}
    λ::U
    T::U
    ω::U
    tmax::U
    Ez::Function
    vanish::Bool
end

(field::LinearField{U})(t::U) where {U<:AbstractFloat} = field.Ez(t)*box(t,field.tmax,field.vanish)
zero(::Type{LinearField{U}}) where {U<:AbstractFloat} = zero(U)

struct TransverseField{U<:AbstractFloat} <: Field{U}
    λ::U
    T::U
    ω::U
    tmax::U
    Ez::Function
    Ex::Function
    vanish::Bool
end

(field::TransverseField{U})(t::U) where {U<:AbstractFloat} = SVector{2,U}(field.Ez(t),field.Ex(t)).*box(t,field.tmax,field.vanish)
zero(::Type{TransverseField{U}}) where {U<:AbstractFloat} = SVector{2,U}(zero(U), zero(U))

function gdd_params(λ_SI::U, τ₀::U, η::U = zero(U);
                    gdd_phase = false) where {U<:AbstractFloat}
    γ = τ₀^2/(8log(2))
    γ² = γ^2
    η² = η^2

    A = √(γ/√(γ² + η²))
    gdd_phase && (A *= exp(im*atan2(-η,γ)/2))

    a = one(U)./4*(γ/(γ² + η²))
    b = one(U)./2*(η/(γ² + η²))
    A,a,b
end

"""`λ_SI` in meters, `I_SI` in W/cm², `tmax` in cycles,
(Fourier-limited) `fwhm` in seconds, `cep` in units of π, `gdd` in
seconds squared. The pulse will be centered at `tmax/2`. If `vanish`
is `true`, the pulse will identically vanish outside the interval
[0,`tmax`]. `ξ` is the ratio of the minor to major axes of the ellipse
(if non-zero, will be 3d pulse automatically, if zero, will be 2d by
default, but can be overridden)."""
function pulse(λ_SI::U, I_SI::U,
               tmax::U, fwhm::U,
               q::U = one(U),
               cep::U = zero(U), gdd::U = zero(U);
               ξ::U = zero(U),
               threeD::Bool = ξ != zero(U),
               gdd_phase = false,
               vanish = true) where {U<:AbstractFloat}
    ξ != 0 && !threeD && error("Can't specify non-linearly polarized pulses in 2d")

    λ,T,ω = fundamental(λ_SI)
    T2 = U(2.41888430e-17T)^2
    I = intensity(I_SI)
    E₀ = √(I)

    A,a,b = gdd_params(λ_SI, fwhm, gdd; gdd_phase = gdd_phase)

    # Transfer to cycles as time base
    a *= T2
    b *= T2

    φ₀ = cep*π

    E = t -> E₀/√(1+ξ^2)*exp(im*2π*q*(t-tmax/2) + im*q*φ₀)*A*exp(-(a-im*b)*(t-tmax/2)^2)

    threeD ?
        TransverseField(λ, T, ω, tmax, t -> U(real(E(t))), t -> U(real(im*ξ*E(t))), vanish) :
        LinearField(λ, T, ω, tmax, t -> U(real(E(t))), vanish)
end

function strong_field_params(λ_SI::U, I_SI::U, Ip::U) where {U<:AbstractFloat}
    λ,T,ω = fundamental(λ_SI)
    I = intensity(I_SI)
    E₀ = √(I)
    Up = I/(4*ω^2)
    cutoff = 3.17Up - Ip
    keldysh = ω*√(2abs(Ip))/E₀
    Dict(:lambda => λ,
         :lambda_SI => λ_SI,
         :T => T,
         :omega => ω,
         :I => I,
         :I_SI => I_SI,
         :E0 => E₀,
         :A0 => E₀/ω,
         :alpha0 => E₀/ω^2,
         :Ip => Ip,
         :Up => Up,
         :keldysh => keldysh,
         :cutoff => cutoff,
         :cutoff_HO => cutoff/ω)
end

export Field, CompositeField, fundamental, delay, gdd_params, pulse, eltype, (+), top_hat, top_hat_trunc, strong_field_params
end # module
