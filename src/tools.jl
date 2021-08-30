"""
Modify pitch using time-scale modification `tsm`.
"""
function pitchshift(tsm::AT, x::AbstractVector{T}, semitones::Real; fs::Real=1) where {T<:Number,AT<:AbstractTimeScaleModifier}
    rate = 2 ^ (semitones / 12)
    y = tsmodify(tsm, x, rate)
    y′ = convert.(T, resample(y, 1 / rate)) # resample returns Float64
    fixlength(y′, length(x))
end

"""
Change speed of `x` without affecting its pitch using time-scale modification `tsm`.
"""
function timestretch(tsm::AT, x::AbstractVector{T}, speed::Real; fs::Real=1) where {T<:Number,AT<:AbstractTimeScaleModifier}
    y = tsmodify(tsm, x, speed)
    fixlength(y, length(x))
end