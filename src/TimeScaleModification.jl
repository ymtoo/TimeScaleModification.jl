module TimeScaleModification

using DSP
using Interpolations
using Reexport

@reexport using DSP.Windows
export OLA, WSOLA, tsmodify
export pitchshift, timestretch

abstract type AbstractTimeScaleModifier end

include("utils.jl")
include("core.jl")
include("tools.jl")

end # module
