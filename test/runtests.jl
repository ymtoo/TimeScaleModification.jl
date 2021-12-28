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
    Ts = [Float64, Float32]
    
    # OLA
    @test OLA(n, synhopsize) == OLA(n, synhopsize, rect)
    ola = OLA(n, synhopsize, hanning)
    @test ola.winfunc == hanning
    @test gettolarance(ola) == 0 
    for T ∈ Ts
        x1 = convert.(T, x)
        for s ∈ [0.5,1.0,1.5]
            @test length(tsmodify(ola, x1, s)) == round(Int, length(x1) * s)
        end
    end

    # WSOLA 
    @test WSOLA(n, synhopsize) == WSOLA(n, synhopsize, rect, 1)
    @test WSOLA(n, synhopsize, hanning) == WSOLA(n, synhopsize, hanning, 1)
    tol = 10
    wsola = WSOLA(n, synhopsize, hanning, tol)
    @test wsola.tolerance == tol
    @test gettolarance(wsola) == convert(eltype(x), tol)
    for T ∈ Ts
        x1 = convert.(T, x)
        for s ∈ [0.5,1.0,1.5]
            y = tsmodify(wsola, x1, s)
            @test length(y) == round(Int, length(x1) * s)
            s == 1.0 && (@test y ≈ x1 atol=1e-3) 
        end
    end

    # PhaseVocoder
    @test PhaseVocoder(n, synhopsize) == PhaseVocoder(n, synhopsize, rect)
    @test PhaseVocoder(n, synhopsize, hanning) == PhaseVocoder(n, synhopsize, hanning, 0, false, false, false)
    phasevocoder = PhaseVocoder(n, synhopsize, hanning)
    for T ∈ Ts
        x1 = convert.(T, x)
        for s ∈ [0.5,1.0,1.5]
            y = tsmodify(phasevocoder, x1, s)
            @test length(y) == round(Int, length(x1) * s)
            s == 1.0 && (@test y ≈ x1 atol=1e-3) 
        end
    end
end

@testset "tools" begin
    
    tsms = [OLA(256, 128, hanning),
           WSOLA(256, 128, hanning, 10), 
           PhaseVocoder(256, 128, hanning, 16, false, false, true)]
    for tsm ∈ tsms
        xps1 = pitchshift(tsm, x, 0)
        xts1 = timestretch(tsm, x, 1)
        @test xps1 ≈ x atol=1e-2
        @test xts1 ≈ x atol=1e-2
    end
end