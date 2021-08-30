"""
Overlap-Add.
"""
struct OLA <: AbstractTimeScaleModifier
    n::Int
    synhopsize::Int
    winfunc::Function
end
OLA(n, synhopsize) = OLA(n, synhopsize, rect)

"""
Waveform Similarity Overlap-Add. `tolerance` is the number of samples the 
window positions in the input signal may be shifted to avoid phase discontinuities
when overlap-adding them to form the output signal.
"""
struct WSOLA <: AbstractTimeScaleModifier
    n::Int
    synhopsize::Int
    winfunc::Function
    tolerance::Int 
end
WSOLA(n, synhopsize, winfunc) = WSOLA(n, synhopsize, winfunc, 1)
WSOLA(n, synhopsize) = WSOLA(n, synhopsize, rect)

gettolarance(::OLA) = 0
gettolarance(f::WSOLA) = f.tolerance

function tsmodify(f::Union{OLA,WSOLA}, x::AbstractVector{T}, s::Real) where {T<:Number}
    tol = gettolarance(f)
    anchorpoints = getanchorpoints(x, s)
    Nout = anchorpoints[end,2]

    win = convert.(T, f.winfunc(f.n))
    winlenhalf = f.n÷2
    synwinpos = 1:f.synhopsize:(Nout + winlenhalf)
    # convert to Float64 due to the overflow of integers for the output ratio 
    anawinpos = LinearInterpolation(Float64.(anchorpoints[:,2]), 
                                    Float64.(anchorpoints[:,1]);
                                    extrapolation_bc=Line())(synwinpos) |>
                x -> round.(Int, convert.(Float64, x)) 
    anahopsize = [0;anawinpos[2:end]-anawinpos[1:end-1]]

    # zero padding
    minfac = minimum(f.synhopsize ./ anahopsize)
    xpad = [zeros(T, winlenhalf + tol);x;zeros(T, ceil(Int, 1/minfac) * f.n + tol)]
    anawinpos .+= tol

    y = zeros(T, Nout + 2 * f.n)
    ow = zeros(T, Nout + 2 * f.n)
    del = 0 # shift of the current analysis window position
    for i = 1:length(anawinpos)-1
        currsynwinran = synwinpos[i]:synwinpos[i]+f.n-1
        curranawinran = (anawinpos[i] + del):(anawinpos[i] + f.n - 1 + del)
        y[currsynwinran] .+= @view(xpad[curranawinran]) .* win
        ow[currsynwinran] .+= win

        natprog = @view(xpad[curranawinran .+ f.synhopsize])

        nextanawinran = (anawinpos[i+1] - tol):(anawinpos[i+1] + f.n - 1 + tol)
        xnext = @view(xpad[nextanawinran])
        cc = fastconv(reverse(xnext), natprog)[f.n:end-f.n+1] # xcorr(xnext, natprog; padmode=:none)[f.n:end-f.n+1]#
        maxindex = argmax(cc)
        del = tol - maxindex + 1
    end
    y[synwinpos[end]:synwinpos[end]+f.n-1] .+= xpad[anawinpos[end]+del:anawinpos[end]+f.n-1+del] .* win
    ow[synwinpos[end]:synwinpos[end]+f.n-1] .+= win
    ow[ow .< 1e-3] .= one(T)
    y ./= ow
    y = y[winlenhalf+1:end]
    y[1:Nout]
end

