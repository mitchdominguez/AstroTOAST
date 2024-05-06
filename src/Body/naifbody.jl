#=
naifbody.jl

Celestial body objects

Adapted from code written by Rolfe Power
=#
"""
    AbstractCelestialBody

Abstract base type for all celestial bodies.
"""
abstract type AbstractCelestialBody end

"""
    NaifBody

Abstract type for for all NAIF body types.
"""
abstract type NaifBody <: AbstractCelestialBody end

"""
    DefaultNaifBody

Concrete type for predefined default NAIF bodies with singleton implementations.
"""
struct DefaultNaifBody <: NaifBody
    name::String
    naif_id::Int
    parent_body::Union{Nothing, NaifBody}
    gravitational_parameter::Float64
    mean_radius::Union{Nothing,Float64}
    semimajor_axis::Union{Nothing,Float64}
    eccentricity::Union{Nothing,Float64}
end

# Define some useful bodies and a dictionary to hold everything in
Bodies = Dict{String,NaifBody}()
Sun = DefaultNaifBody("Sun", # Name
                        10,  # NAIF_ID
                        nothing, # Parent Body
                        1.3271244004193930e+11, # Gravitational Parameter
                        6.9600000000000000e+05, # Mean Radius
                        nothing,  # Semimajor Axis
                        nothing) # Eccentricity
merge!(Bodies, Dict("Sun"=>Sun))


Earth = DefaultNaifBody("Earth", # Name
                        399,  # NAIF_ID
                        Sun, # Parent Body
                        3.9860043543609593e+05, # Gravitational Parameter
                        6.3710083666666660e+03, # Mean Radius
                        1.4959789217545033e+08,  # Semimajor Axis
                        1.6735932113458880e-02) # Eccentricity
merge!(Bodies, Dict("Earth"=>Earth))

Moon = DefaultNaifBody("Moon", # Name
                        302,  # NAIF_ID
                        Earth, # Parent Body
                        4.9028000661637961e+03, # Gravitational Parameter
                        1.7374000000000003e+03, # Mean Radius
                        3.8474799201129237e+05,  # Semimajor Axis
                        9.0000558455574642e-02) # Eccentricity
merge!(Bodies, Dict("Moon"=>Moon))

EMB = DefaultNaifBody("Earth-Moon Barycenter", # Name
                        3,  # NAIF_ID
                        Sun, # Parent Body
                        4.0350323550225981e+05, # Gravitational Parameter
                        nothing,
                        1.4959789401764473e+08,  # Semimajor Axis
                        1.6729203447508743e-02) # Eccentricity
merge!(Bodies, Dict("Earth-Moon Barycenter"=>EMB))

"""
    Base.show

Overload the show operator to pretty print the NAIF body to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", cb::DefaultNaifBody)
    print(io, """Celestial Body $(cb.name):
    \tNAIF ID: $(cb.naif_id)
    \tGM     : $(cb.gravitational_parameter) km^3/s^2""")
    if typeof(cb.mean_radius) <: Float64
        print(io, """ \n\tRadius : $(cb.mean_radius) km""")
    end
    if typeof(cb.parent_body) <: NaifBody
        print(io, "\n\tParent : $(cb.parent_body.name)" )
    end
end
