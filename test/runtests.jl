using TimeScaleModification
using TimeScaleModification: getanchorpoints, fixlength, gettolarance, fastconv, _stft, _istft 

using DSP
using SignalAnalysis
using Test

fs = 96000
freq = 5000
x = cw(freq, 0.5, fs) |> real |> collect
N = length(x)

@testset "utils" begin

    s1 = 0.9
    @test getanchorpoints(x, s1) == [1 1; N ceil(Int, s1 * N)]
    s2 = [1 1; N N÷2]
    @test getanchorpoints(x, s2) == s2

    y1 = fixlength(x, N÷2)
    @test length(y1) == N÷2
    @test y1 == x[1:N÷2]
    y2 = fixlength(x, 2N)
    @test length(y2) == 2N
    @test y2[1:N] == x
    @test y2[N+1:end] == zeros(eltype(x),N)

    a = randn(256)
    b = randn(256)
    @test conv(a, b) ≈ fastconv(a, b) atol=1e-9

    n = 256
    synhopsize = 128
    spec, f, t = _stft(x, n, synhopsize; fs=fs)
    freqindex = argmin(abs.(f .- freq))
    for xcol ∈ eachcol(spec)
        if !all(xcol .≈ 0)
            @test argmax(abs.(xcol)) == freqindex
        end
    end

    for isfftshift ∈ [true,false]
        for zeropad ∈ 0:10
            spec, f, t = _stft(x, n, synhopsize; fs=fs, zeropad=zeropad, isfftshift=isfftshift)
            for numiter ∈ 1:10
                for isrestore ∈ [true,false]
                    xt1 = _istft(spec, n, synhopsize; zeropad=zeropad, numiter=numiter, origsiglen=length(x), isfftshift=isfftshift, isrestore=isrestore)
                    win = rect(n)
                    xt2 = _istft(spec, win, synhopsize; zeropad=zeropad, numiter=numiter, origsiglen=length(x), isfftshift=isfftshift, isrestore=isrestore)
                    @test xt1 == xt2 
                    @test xt1 ≈ x atol=1e-9
                end
            end
        end
    end
end

@testset "core" begin
    
    n = 256
    synhopsize = 128
    
    # OLA
    @test OLA(n, synhopsize) == OLA(n, synhopsize, rect)
    ola = OLA(n, synhopsize, hanning)
    @test ola.winfunc == hanning
    @test gettolarance(ola) == 0 
    for s ∈ [0.5,1.0,1.5]
        @test length(tsmodify(ola, x, s)) == round(Int, length(x) * s)
    end

    # WSOLA 
    @test WSOLA(n, synhopsize) == WSOLA(n, synhopsize, rect, 1)
    @test WSOLA(n, synhopsize, hanning) == WSOLA(n, synhopsize, hanning, 1)
    tol = 10
    wsola = WSOLA(n, synhopsize, hanning, tol)
    @test wsola.tolerance == tol
    @test gettolarance(wsola) == convert(eltype(x), tol)
    for s ∈ [0.5,1.0,1.5]
        y = tsmodify(wsola, x, s)
        @test length(y) == round(Int, length(x) * s)
        s == 1.0 && (@test y ≈ x atol=1e-9) 
    end

    # PhaseVocoder
    @test PhaseVocoder(n, synhopsize) == PhaseVocoder(n, synhopsize, rect)
    @test PhaseVocoder(n, synhopsize, hanning) == PhaseVocoder(n, synhopsize, hanning, 0, false, false, false)
    phasevocoder = PhaseVocoder(n, synhopsize, hanning)
    for s ∈ [0.5,1.0,1.5]
        y = tsmodify(phasevocoder, x, s)
        @test length(y) == round(Int, length(x) * s)
        s == 1.0 && (@test y ≈ x atol=1e-9) 
    end

end

@testset "tools" begin
    


end