# Functions dealing with storing and recovering AstroTOAST data 
# structures using external files 
using AstroTOAST
using MAT


# Abstract types for ease of referring to AstroTOAST types
abstract type ATType end

struct _PeriodicOrbit <: ATType end


function from_dict(dd::Dict)
    # Establish data type of data saved in the Dict
    datatype = get(dd, "datatype", nothing)

    if isnothing(datatype)
        throw(ErrorException("No data type is specified in the input dictionary!"))
    end

    ### TODO check that the rest of the entries can be used to instantiate the struct
    if datatype=="PeriodicOrbit"
        return from_dict(dd, eval(Symbol("_"*datatype))())
    end

    # return from_dict(dd::Dict, datatype)
end

function from_dict(dd::Dict, ::_PeriodicOrbit)
    ## Extract model
    _P1 = get(dd, "P1", nothing)
    _P2 = get(dd, "P2", nothing)
    _mu = get(dd, "mu", nothing)

    # Check that data is input correctly
    if isnothing(_P1) || isnothing(_P2) || isnothing(_mu)
        throw(ErrorException("One of P1, P2, or mu does not exist! A dynamical model cannot be specified"))
    end

    # Create model
    if isempty(_P1)
        model = Cr3bpModel(_mu)
    else
        model = Cr3bpModel(Bodies[_P1], Bodies[_P2])
    end

    ## Extract PO parameters
    _name = get(dd, "name", nothing)
    _family = get(dd, "family", nothing)
    _X0 = get(dd, "X0", nothing)
    _TOF = get(dd, "TOF", nothing)
    if isnothing(_name) || isnothing(_family) || isnothing(_X0) || isnothing(_TOF)
        throw(ErrorException("One of name, family, or X0, or TOF does not exist! A dynamical model cannot be specified"))
    end

    # Try to create the periodic orbit as is from the dictionary
    # If that fails, then retarget using continuity
    try 
        return PeriodicOrbit(model, _X0, _TOF, _name, _family)
    catch
        po = targetcontinuity(model, _X0, _TOF, _name, _family)[1]
        println("Converged")
        return po
    end
end

"""
    from_mat(matfile::String)

Return objects stored in `matfile`. The assumed format is that `matread(matfile)` returns a
Dict which at the top level has an entry `data`. Within this entry is another Dict

"""
function from_mat(matfile::String)
    matdata = matread(matfile)

    # show(stdout, "text/plain", matdata)
    ### TODO I need to specify a general data format for everything saved to a mat file from AstroTOAST 
    
    datadictexists = ~isnothing(get(matdata, "data", nothing))

    outobj = []

    for i = 1:length(matdata["data"])
        push!(outobj, from_dict(matdata["data"][i]))
    end

    return outobj
end
