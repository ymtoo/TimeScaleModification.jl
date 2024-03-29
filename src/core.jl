"""
Overlap-Add. `n` is window size, `synhopsize` is the hop size of the synthesis window and 
`winfunc` is the analysis and synthesis window for STFT.
"""
struct OLA{W} <: AbstractTimeScaleModifier
    n::Int
    synhopsize::Int
    winfunc::W
end
OLA(n, synhopsize) = OLA(n, synhopsize, rect)

"""
Waveform Similarity Overlap-Add. `n` is window size, `synhopsize` is the hop size of the synthesis window and 
`winfunc` is the analysis and synthesis window for STFT.

`tolerance` is the number of samples the window positions in the input signal may be shifted 
to avoid phase discontinuities when overlap-adding them to form the output signal.
"""
struct WSOLA{W} <: AbstractTimeScaleModifier
    n::Int
    synhopsize::Int
    winfunc::W
    tolerance::Int 
end
WSOLA(n, synhopsize, winfunc) = WSOLA(n, synhopsize, winfunc, 1)
WSOLA(n, synhopsize) = WSOLA(n, synhopsize, rect)

gettolarance(::OLA) = zero(Int)
gettolarance(f::WSOLA) = f.tolerance

"""
Modify time scale of `x` using either `OLA` or `WSOLA`.

# Arguments
- tsm: OLA or WSOLA instance
- x: 1-D signal
- s: a constant scaling factor or a n x 2 matrix representing a set of n anchorpoints 
relating sample positions in the input signal with sample positions in the output signal.

# Returns
Time scale modifed signal

# Examples
```julia-repl
julia> tsmodify(OLA(256,128,hanning), randn(1000), 1.5)
1500-element Vector{Float64}:
 -0.5010944306688007
  2.0751141429688493
  ⋮
  0.8717305567708631
 -0.9940213359357077

julia> tsmodify(WSOLA(256,128,hanning,10), randn(1000), [1 1; 10 10])
10-element Vector{Float64}:
  0.07891842073223751
 -0.11009242947582566
  ⋮
 -2.342067676109941
  0.2866198482301998
```
"""
function tsmodify(tsm::Union{OLA,WSOLA}, x::AbstractVector{T}, s::Union{Real,AbstractMatrix{Int}}) where {T<:Number}
    anchorpoints = getanchorpoints(x, s)
    Nout = anchorpoints[end,2]

    win = convert.(T, tsm.winfunc(tsm.n))
    winlen = length(win)
    winlenhalf = winlen ÷ 2

    synwinpos = 1:tsm.synhopsize:(Nout + winlenhalf)
    # convert to Float64 due to the overflow of integers for the output ratio 
    anawinpos = linear_interpolation(Float64.(anchorpoints[:,2]), 
                                     Float64.(anchorpoints[:,1]);
                                     extrapolation_bc=Line())(synwinpos) |>
                x -> round.(Int, convert.(Float64, x)) 
    anahopsize = [0;anawinpos[2:end]-anawinpos[1:end-1]]

    tol = gettolarance(tsm)

    # zero padding
    minfac = minimum(tsm.synhopsize ./ anahopsize)
    xpad = [zeros(T, winlenhalf + tol);x;zeros(T, ceil(Int, 1/minfac) * tsm.n + tol)]
    anawinpos .+= tol

    y = zeros(T, Nout + 2 * winlen)
    ow = zeros(T, Nout + 2 * winlen)
    del = zero(Int) # shift of the current analysis window position
    for i = 1:length(anawinpos)-1
        currsynwinran = synwinpos[i]:(synwinpos[i] + winlen - 1)
        curranawinran = (anawinpos[i] + del):(anawinpos[i] + winlen - 1 + del)
        @views y[currsynwinran] .+= xpad[curranawinran] .* win
        @views ow[currsynwinran] .+= win

        @views natprog = xpad[curranawinran .+ tsm.synhopsize]

        nextanawinran = (anawinpos[i+1] + winlen - 1 + tol):-1:(anawinpos[i+1] - tol)
        @views xnext = xpad[nextanawinran]
        @views cc = fastconv(xnext, natprog)[winlen:end-winlen+1] # xcorr(xnext, natprog; padmode=:none)[f.n:end-f.n+1]#
        del = tol - argmax(cc) + 1
    end
    @views y[synwinpos[end]:synwinpos[end]+winlen-1] .+= xpad[anawinpos[end]+del:anawinpos[end]+tsm.n-1+del] .* win
    @views ow[synwinpos[end]:synwinpos[end]+winlen-1] .+= win
    ow[ow .< T(1e-3)] .= one(T)
    y ./= ow
    y = y[winlenhalf+1:end]
    y[1:Nout]
end

"""
Phase vocoder. `n` is window size, `synhopsize` is the hop size of the synthesis window and 
`winfunc` is the analysis and synthesis window for STFT.

`zeropad` is the number of zeros to be padded to the window. If `isfftshift` is `true`, the 
zero-frequency component is shifted to the center. If `isrestore` is `true`, every windowed 
synthesis frame is rescaled to compensate for energy leakage. If `isphaselock` is `true`, 
phase locking is applied. Details can be found at 
https://www.audiolabs-erlangen.de/content/resources/MIR/TSMtoolbox/pvTSM.m
"""
struct PhaseVocoder{W} <: AbstractTimeScaleModifier
    n::Int
    synhopsize::Int
    winfunc::W
    zeropad::Int
    isfftshift::Bool
    isrestore::Bool
    isphaselock::Bool
end
PhaseVocoder(n, synhopsize, winfunc) = PhaseVocoder(n, synhopsize, winfunc, 0, false, false, false)
PhaseVocoder(n, synhopsize) = PhaseVocoder(n, synhopsize, rect)

"""
Modify time scale of `x` using `PhaseVocoder`.

# Arguments
- tsm: PhaseVocoder instance
- x: 1-D signal
- s: a constant scaling factor or a n x 2 matrix representing a set of n anchorpoints 
relating sample positions in the input signal with sample positions in the output signal.

# Returns
Time scale modifed signal

# Examples:
```julia-repl
julia> tsmodify(PhaseVocoder(256,128,hanning,16,false,false,true), randn(1000), 1.5)
1500-element Vector{Float64}:
  0.6042886457594456
 -0.4197586697067266
  ⋮
  0.15547596878459596
 -0.39570063285685614

julia> tsmodify(PhaseVocoder(256,128,hanning,16,false,false,true), randn(1000), [1 1; 10 10])
10-element Vector{Float64}:
  0.9705544964709033
 -0.42565730762449744
  ⋮
  0.7656078049589671
 -0.354979106782902
```
"""
function tsmodify(tsm::PhaseVocoder, 
                  x::AbstractVector{T}, 
                  s::Union{Real,AbstractMatrix{Int}}) where {T<:Number}
    anchorpoints = getanchorpoints(x, s)
    Nout = anchorpoints[end,2]

    win = convert.(T, tsm.winfunc(tsm.n)) 
    #win = [zeros(T, tsm.zeropad÷2);x;zeros(T, tsm.zeropad÷2)]
    #winlen = length(win)
    #winlenhalf = winlen ÷ 2
    winlen = tsm.n + tsm.zeropad
    winlenhalf = winlen ÷ 2

    synwinpos = 1:tsm.synhopsize:(Nout + winlenhalf)
    # convert to Float64 due to the overflow of integers for the output ratio 
    anawinpos = linear_interpolation(Float64.(anchorpoints[:,2]), 
                                     Float64.(anchorpoints[:,1]);
                                     extrapolation_bc=Line())(synwinpos) |>
                x -> round.(Int, convert.(Float64, x)) 
    anahopsize = [0;anawinpos[2:end]-anawinpos[1:end-1]]

    #y = zeros(T, Nout)
    spec, f, t = _stft(x, win, anawinpos; zeropad=tsm.zeropad, isfftshift=tsm.isfftshift)
    Y = zero(spec)
    Y[:,1] = spec[:,1]
    
    k = 0:winlenhalf
    ω = 2π .* k ./ winlen
    
    for i ∈ axes(spec, 2)[2:end] #2:size(spec,2)
        δϕ = ω * anahopsize[i]
        ϕcurr = angle.(spec[:,i])
        ϕlast = angle.(spec[:,i-1])
        hpi = (ϕcurr - ϕlast) .- δϕ
        hpi .-= 2π .* round.(hpi ./ (2π)) 
        
        ipasample = ω .+ hpi ./ anahopsize[i]
        ipahop = ipasample .* tsm.synhopsize

        ϕsyn = angle.(Y[:,i-1])

        if !tsm.isphaselock
            phasor = exp.(im .* (ϕsyn .+ ipahop .- ϕcurr))
        else
            pkindices, irs, ire = findpeaks(abs.(spec[:,i]))
            θ = zero(Y[:,i])
            for n ∈ eachindex(pkindices) #1:length(pkindices)
                @views θ[irs[n]:ire[n]] .= ϕsyn[pkindices[n]] .+ ipahop[pkindices[n]] .- ϕcurr[pkindices[n]]
            end
            phasor = exp.(im .* θ)
        end
        @views Y[:,i] = phasor .* spec[:,i]
    end
    _istft(Y, 
           win, 
           tsm.synhopsize; 
           zeropad=tsm.zeropad, 
           numiter=1, 
           origsiglen=Nout,
           isfftshift=tsm.isfftshift,
           isrestore=tsm.isrestore)
end
