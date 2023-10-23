using MATLAB

export get_odesolver, get_rhsfunction, get_abstol, get_reltol

struct MatlabODESolver
    odesolver::String
    rhsfunction::String
    abstol::Union{Float64, Vector{Float64}}
    reltol::Union{Float64, Vector{Float64}}

    function MatlabODESolver(odesolver::String, rhsfunction::String; abstol=[], reltol=[])
        if rhsfunction[1] == '@'
            rhsstr = rhsfunction
        else
            rhsstr = "@"*rhsstr
        end

        return new(odesolver, rhsstr, abstol, reltol)
    end
end

get_odesolver(mos::MatlabODESolver) = mos.odesolver
get_rhsfunction(mos::MatlabODESolver) = mos.rhsfunction
get_abstol(mos::MatlabODESolver) = mos.abstol
get_reltol(mos::MatlabODESolver) = mos.reltol

function matlabsolve(mos::MatlabODESolver, q0::Vector{Float64}, times::AbstractVector)
    str = "[t,Z] = $(get_odesolver(mos))($(get_rhsfunction(mos)), $(times), $(q0), odeset('AbsTol', $(get_abstol(mos)), 'RelTol', $(get_reltol(mos))));"
    eval_string(str)
    @mget t
    @mget Z
    return t, Z
end

