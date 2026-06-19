using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using Plots
gr()

include(joinpath(@__DIR__, "..", "src", "RTNHE.jl"))
using .RTNHE

function main(args = ARGS)
    if isempty(args)
        println("No mode provided. Defaulting to identical.")
        RTNHE.main(["identical"])
    else
        RTNHE.main(args)
    end
end

main(ARGS)