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

"""
Computes the STFT of a signal `x`.
"""
function _stft(x::AbstractVector{T}, 
               win::AbstractVector{T},
               anahopsize::Union{Int,Vector{Int}};
               fs::Real=1,
               zeropad::Int=0,
               numframes::Union{Nothing,Int}=nothing,    
               isfftshift::Bool=false) where {T<:Number}
    win = [zeros(T, zeropad÷2);win;zeros(T, zeropad÷2)]
    winlen = length(win)
    winlenhalf = winlen÷2

    # zero padding
    maxanahopsize = maximum(anahopsize)
    xpad = [zeros(T, winlenhalf);x;zeros(T, winlen + maxanahopsize)]

    if length(anahopsize) == 1
        numframes === nothing && (numframes = floor(Int, (length(xpad) - winlen) / anahopsize + 1))
        winpos = (0:numframes-1) .* anahopsize .+ 1
    else
        numframes === nothing && (numframes = length(anahopsize))
        winpos = anahopsize[1:numframes]
    end

    spec = T <: Real ? zeros(Complex{T}, winlenhalf+1, numframes) : zeros(T, winlenhalf+1, numframes)
    for i ∈ 1:numframes
        xi = xpad[winpos[i]:winpos[i] + winlen - 1] .* win
        isfftshift && (xi = DSP.fftshift(xi))
        Xi = DSP.fft(xi)
        spec[:,i] = Xi[1:winlenhalf+1]
    end
    t = (winpos .- 1) ./ fs
    f = (0:winlenhalf) .* fs ./ winlen
    spec, f, t
end
function _stft(x::AbstractVector{T}, 
               n::Int, 
               anahopsize::Union{Int,Vector{Int}};
               fs::Real=1,
               winfunc::Function=rect,
               zeropad::Int=0,
               numframes::Union{Nothing,Int}=nothing,  
               isfftshift::Bool=false) where {T<:Number}
    win = convert.(T, winfunc(n))
    _stft(x, win, anahopsize; fs=fs, zeropad=zeropad, numframes=numframes, isfftshift=isfftshift)
end

"""
Computes the ISTFT of `spec`.

# Reference
D. W. Griffin and J. S. Lim, "Signal estimation from modified short-time Fourier transform," IEEE Trans. ASSP, vol.32, no.2, pp.236–243, Apr. 1984 
"""
function _istft(spec::AbstractMatrix{Complex{T}}, 
                n::Int, 
                synhopsize::Int;
                winfunc::Function=rect,
                zeropad::Int=0,
                numiter::Int=1,
                origsiglen::Union{Nothing,Int}=nothing,
                isfftshift::Bool=false,
                isrestore::Bool=false) where {T<:Real}
    numframes = size(spec, 2)

    # first iteration
    yi = lsee_mstft(spec, n, synhopsize; winfunc=winfunc, zeropad=zeropad, isfftshift=isfftshift, isrestore=isrestore)
    for j ∈ 2:numiter
        Yi = abs.(spec) .* exp.(im .* angle.(first(_stft(yi, n, synhopsize; winfunc=winfunc, zeropad=zeropad, numframes=numframes, isfftshift=isfftshift))))
        yi = lsee_mstft(Yi, n, synhopsize; winfunc=winfunc, zeropad=zeropad, isfftshift=isfftshift, isrestore=isrestore)
    end
    if origsiglen === nothing
        yi
    else
        yi[1:origsiglen]
    end
end
function lsee_mstft(spec::AbstractMatrix{Complex{T}}, 
                    n::Int, 
                    synhopsize::Int;
                    #fs::Real=1,
                    winfunc::Function=rect,
                    zeropad::Int=0,
                    isfftshift::Bool=false,
                    isrestore::Bool=false
                    ) where {T<:Real}
    win = convert.(T, winfunc(n))
    win = [zeros(T, zeropad÷2);win;zeros(T, zeropad÷2)]
    winlen = length(win)
    winlenhalf = winlen÷2

    numframes = size(spec, 2)
    winpos = (0:numframes-1) .* synhopsize .+ 1
    xlen = winpos[end] + winlen - 1

    x = zeros(T, xlen)
    ow = zeros(T, xlen)
    for i ∈ 1:numframes
        currspec = spec[:,i]
        Xi = [currspec;conj.(currspec[end-1:-1:2])]
        xi = real.(DSP.ifft(Xi))
        isfftshift && (xi = DSP.fftshift(xi))
        xiw = xi .* win

        if isrestore
            xienergy = sum(abs, xi)
            xiwenergy = sum(abs, xiw)
            xiw .*= xienergy ./ (xiwenergy .+ eps(T))
        end

        x[winpos[i]:winpos[i]+winlen-1] .+= xiw
        ow[winpos[i]:winpos[i]+winlen-1] .+= win.^2  
    end
    ow[ow .< 1e-3] .= one(T)
    x ./= ow
    x[winlenhalf+1:end-winlenhalf]
end