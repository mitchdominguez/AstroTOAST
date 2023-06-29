using AstroTOAST

abstract type AbstractOrbitalElements end

"""
    struct ClassicalOrbitalElements

Struct containing the Keplerian (classical) orbital elements: semimajor axis,
eccentricity, right ascension of the ascending node, inclination, and argument
of periapsis. All angular orbital elements are stored internally as radians.
"""
struct ClassicalOrbitalElements
    coe::Vector{Float64} 

    """
        function ClassicalOrbitalElements(a::Real, e::Real, i::Real, Ω::Real, ω::Real; degrees=true)   

        Constructor for ClassicalOrbitalElements. `a` in km, `e` nondimensional, and 
        `i`, `Ω`, and `ω` are assumed to be degrees by default. Set `degrees=false` to
        input radians.
    """
    function ClassicalOrbitalElements(a::Real, e::Real, Ω::Real, i::Real, ω::Real; degrees=true)
        if degrees
            # Angular orbital elements given in degrees

            # Check that inclination is legal
            if i>90 || i<-90
                DomainError("Inclination must be between -90 and 90 degrees") |> throw
            end
            
            new([a,e,(deg2rad.([Ω,i,ω]).%(2pi))...])
        else

            # Check that inclination is legal
            if i>pi || i<-pi
                DomainError("Inclination must be between -π and π radians") |> throw
            end
            
            # Angular orbital elements given in radians
            new([a,e,Ω%(2pi),i,ω%(2pi)])
        end
    end
end

"""
    sma(coe::ClassicalOrbitalElements)

Return the semimajor axis
"""
sma(coe::ClassicalOrbitalElements) = coe.coe[1]

"""
    ecc(coe::ClassicalOrbitalElements)

Return the eccentricity
"""
ecc(coe::ClassicalOrbitalElements) = coe.coe[2]

"""
    raan(coe::ClassicalOrbitalElements)

Return the RAAN. Default return value in radians. Specify `degrees=true` for degrees
"""
raan(coe::ClassicalOrbitalElements; degrees=false) = degrees ? rad2deg(coe.coe[3]) : coe.coe[3]

"""
    inc(coe::ClassicalOrbitalElements)

Return the inclination. Default return value in radians. Specify `degrees=true` for degrees
"""
inc(coe::ClassicalOrbitalElements ; degrees=false) = degrees ? rad2deg(coe.coe[4]) : coe.coe[4]

"""
    aop(coe::ClassicalOrbitalElements)

Return the argument of periapsis. Default return value in radians. Specify `degrees=true` for degrees
"""
aop(coe::ClassicalOrbitalElements ; degrees=false) = degrees ? rad2deg(coe.coe[5]) : coe.coe[5]
