"""
Modify pitch using time-scale modifier `tsm`.

# Arguments
- tsm: time-scale modifier 
- x: 1-D signal
- semitones: how many semitones to shift
- fixlen: if `true`, right padding the output with zeros such that the lengths of the output and `x` are the same. 

# Returns
Pitch-shifted signal

# Examples
```julia-repl
julia> pitchshift(OLA(256,128,hanning), randn(1000), 3; fixlen=false)
1001-element Vector{Float64}:
 -0.4600176073023945
  0.158860158951877
  ⋮
  0.22459122115497698
  0.33112956560849627

julia> pitchshift(OLA(256,128,hanning), randn(1000), -3; fixlen=true)
1000-element Vector{Float64}:
  -0.3279038700745687
  0.07825216616901506
  ⋮
  -1.3089405766166855
  -0.21792912756346938
```
"""
function pitchshift(tsm::AT, 
                    x::AbstractVector{T}, 
                    semitones::Real; 
                    fixlen::Bool=true) where {T<:Number,AT<:AbstractTimeScaleModifier}
    rate = 2 ^ (semitones / 12)
    y = tsmodify(tsm, x, rate)
    y′ = convert.(T, resample(y, 1 / rate)) # resample returns Float64
    fixlen ? fixlength(y′, length(x)) : y′
end

"""
Change speed of `x` without affecting its pitch using time-scale modifier `tsm`.

# Arguments
- tsm: time-scale modifier 
- x: 1-D signal
- speed: if greater than 1, the signal is sped up and if less than 1, the signal is slowed down
- fixlen: if `true`, right padding the output with zeros such that the lengths of the output and `x` are the same. 

# Returns
Time-stretched signal

# Examples
```julia-repl
julia> timestretch(OLA(256,128,hanning), randn(1000), 3; fixlen=false)
3000-element Vector{Float64}:
  0.7586887836708426
  0.7820819246695768
  ⋮
 -0.3460853693317422
 -0.5682688881480474

julia> timestretch(OLA(256,128,hanning), randn(1000), 0.5; fixlen=true)
1000-element Vector{Float64}:
  1.0081982073491438
  1.165989335000566
  ⋮
  0.0
  0.0
```
"""
function timestretch(tsm::AT, 
                     x::AbstractVector{T}, 
                     speed::Real; 
                     fixlen::Bool=true) where {T<:Number,AT<:AbstractTimeScaleModifier}
    y = tsmodify(tsm, x, speed)
    fixlen ? fixlength(y, length(x)) : y
end
