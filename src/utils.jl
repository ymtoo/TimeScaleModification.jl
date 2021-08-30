"""
Convert the scaling factor to anchor points. An anchor point is a pair of sample positions where the first
entry specifies a sample position in the input signal and the second entry specifies a sample position in
the output signal.
"""
function getanchorpoints(x::AbstractVector{T}, s::Union{Real,AbstractMatrix{Int}}) where {T<:Number}
    if length(s) == 1
        return [1 1; size(x,1) ceil(Int, s * size(x,1))]
    elseif (size(s,2) == 2) && (ndims(s) == 2)
        return s
    else
        throw(ArgumentError("Invalid time stretching factor, 
                             use scalar or `n` pairs of input output sample points with dimensions of (n, 2)."))
    end
end

"""
Fix the length `x` to `n`. If `length(x) > n`, `x` is padded with trailing zeros.  
"""
fixlength(x::AbstractVector{T}, n::Int) where {T<:Number} = length(x) > n ? x[1:n] : [x;zeros(T, n-length(x))]

"""
Current implementation of DSP.conv is slow. See https://github.com/JuliaDSP/DSP.jl/issues/167.
"""
function fastconv(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Number}
    m = length(x)
    n = length(y)
    xc = zeros(T, m+n-1)
    @inbounds for j=1:m
        for k=1:n
            xc[j+k-1] += x[j] * y[k]
        end
    end
    xc
end