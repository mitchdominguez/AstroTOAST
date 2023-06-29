using AstroTOAST


#### EXPORT
export TrueAnomaly, MeanAnomaly, EccentricAnomaly

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

###################################
# RETURN CONSTANT ORBITAL ELEMENTS
###################################

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


########################################
# RETURN TIME-VARYING ORBITAL ELEMENTS
########################################

# Type that captures the different kinds of anomalies (true, mean, eccentric)
abstract type Anomaly end;

struct TrueAnomaly
    ta::Float64

    #### Basic Constructor
    function TrueAnomaly(ta::Real; degrees=true)
        if degrees
            # Input in degrees
            new(deg2rad(ta)%(2pi))
        else
            # Input in radians
            new(ta%(2pi))
        end
    end
end

struct MeanAnomaly
    M::Float64
    
    #### Basic Constructor
    function MeanAnomaly(M::Real; degrees=true)
        if degrees
            # Input in degrees
            new(deg2rad(M)%(2pi))
        else
            # Input in radians
            new(M%(2pi))
        end
    end

end

struct EccentricAnomaly
    E::Float64
    
    function EccentricAnomaly(E::Real; degrees=true)
        if degrees
            # Input in degrees
            new(deg2rad(E)%(2pi))
        else
            # Input in radians
            new(E%(2pi))
        end
    end
end

value(TA::TrueAnomaly) = TA.ta
value(M::MeanAnomaly) = M.M
value(E::EccentricAnomaly) = E.E


#### Mean anomaly from eccentric anomaly, orbital elements, i.e. Kepler's eqn
"""
    function MeanAnomaly(E::EccentricAnomaly, coe::ClassicalOrbitalElements)

Return the mean anomaly, given eccentric anomaly and orbital elements, i.e. 
return the mean anomaly from Kepler's equation [M = E-e*sin(E)]
"""
function MeanAnomaly(E::EccentricAnomaly, coe::ClassicalOrbitalElements)
    MeanAnomaly(value(E)-ecc(coe)*sin(value(E)), degrees=false)
end

###################################
# CONVERT BETWEEN COE AND STATES
###################################
