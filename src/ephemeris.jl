## Set spice kernels to use and furnish them
defaultLSK = "data/naif0012.tls"
defaultSPK = "data/de440.bsp"

function load_default_kernels()
    furnsh(defaultLSK)
    furnsh(defaultSPK)
end

function currently_loaded_kernels(;ktype="ALL")
    N = ktotal(ktype)
    if N == 0
        println("No kernels of type '$(ktype)' are currently loaded")
    else
        println("$(N) kernels of type '$(ktype)' are currently loaded")
        println("---")

        for i = 1:N
            fname, ftype = kdata(i, ktype)[1:2]
            println(i, " - ", fname, " - ",  ftype)
        end
    end
end

function clear_all_kernels()
    N = ktotal("ALL")
    for i = 1:N
        fname = kdata(1, "ALL")[1]
        unload(fname)
    end
end

function load_only_default_kernels()
    clear_all_kernels()
    load_default_kernels()
end
