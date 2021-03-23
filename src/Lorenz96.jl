module Lorenz96

    export L96

    include("rhs.jl")
    include("time_integration.jl")
    include("L96.jl")

end
