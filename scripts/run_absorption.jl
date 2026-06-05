using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using Plots
gr()

include(joinpath(@__DIR__, "..", "src", "RTNAbsorption.jl"))
using .RTNAbsorption

function main(args = ARGS)
    if isempty(args)
        println("No mode provided. Defaulting to identical.")
        RTNAbsorption.main(["identical"])
    else
        RTNAbsorption.main(args)
    end
end

main(ARGS)