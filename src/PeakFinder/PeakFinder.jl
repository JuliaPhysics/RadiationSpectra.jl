function peakfinder(source::AbstractArray; σ::Real=2., threshold::Real=10., backgroundRemove::Bool=true, deconIterations::Int=3, markov::Bool=true, averWindow::Int=3)
    sigma::Float64 = σ
    PEAK_WINDOW::Int = 1024 # const of TSpectrum
    if (sigma < 1) throw(ArgumentError("Invalid sigma, must be greater than or equal to 1")) end
    if (threshold<=0 || threshold>=100) throw(ArgumentError("Invalid threshold, must be between 0 and 100")) end
    if markov if (averWindow <= 0) throw(ArgumentError("Averaging Window must be positive")) end end
    ssize = length(source)
    numberIterations = round(Int, 7*sigma + 0.5, RoundDown)
    if backgroundRemove if (ssize < 2*numberIterations + 1) throw(ArgumentError("Too large clipping window")) end end
    j = round(Int, 5*sigma + 0.5, RoundDown)
    if (j >= PEAK_WINDOW/2) throw(ArgumentError("Too large sigma")) end
    peak_index::Int=0; size_ext::Int=ssize+2*numberIterations; shift::Int=numberIterations; bw::Int=2
    plocha::Float64 = 0
    m0low, m1low, m2low, l0low, l1low = zeros(Float64, 5)

    k::Int = round(Int, 2*sigma+0.5, RoundDown)
    if k >= 2
        for i in 0:k-1
            a = i; b=source[i+1]
            m0low += 1; m1low += a; m2low += a^2; l0low += b; l1low += a*b
        end
        detlow = m0low*m2low - m1low*m1low
        l1low = (detlow != 0) ? (-l0low*m1low + l1low*m0low)/detlow : l1low = 0.
        if (l1low > 0) l1low=0. end
    else
        l1low = 0.
    end

    i::Int = 2*round(Int, 7*sigma+0.5, RoundDown)
    working_space::Array{Float64, 1} = zeros(Float64, 7*(ssize+i))
    for i in 0:size_ext-1
        if i < shift
            a = i - shift
            working_space[i+size_ext+1] = source[1]+l1low*a
            if (working_space[i+size_ext+1]<0) working_space[i+size_ext+1] = 0. end
        elseif i >= ssize+shift
            a = i - (ssize - 1 + shift)
            working_space[i+size_ext+1] = source[ssize]
            if (working_space[i+size_ext+1]<0)  working_space[i+size_ext+1] = 0. end
        else
            working_space[i+size_ext+1] = source[i-shift+1]
        end
    end


    if backgroundRemove
        for i in 1:numberIterations
            for j in i:size_ext-i-1
                if !markov
                    a = Float64(working_space[size_ext+j+1])
                    b = (working_space[size_ext+j-i+1] + working_space[size_ext+j+i+1]) / 2.0
                    if (b < a) a=b end
                    working_space[j+1]=a
                else
                    a = working_space[size_ext+j+1]
                    av = 0.
                    men = 0.
                    for w = j-bw:j+bw
                        if ((w>=0) && (w<size_ext))
                            av += working_space[size_ext+w+1]
                            men += 1
                        end
                    end
                    av = av/men
                    b = 0.
                    men = 0.
                    for w in j-i-bw:j-i+bw
                        if ((w>=0) && (w<size_ext))
                             b += working_space[size_ext+w+1]
                             men +=1
                        end
                    end
                    b = b/men
                    c = 0.
                    men = 0.
                    for w in j+i-bw:j+i+bw
                        if ((w>=0) && (w<size_ext))
                             c += working_space[size_ext+w+1]
                             men += 1
                        end
                    end
                    c = c/men
                    b = (b+c)/2.
                    if (b<a) av=b end
                    working_space[j+1]=av;
                end
            end
            for j in i:size_ext-i-1
                working_space[size_ext + j + 1] = working_space[j+1];
            end
        end

        for j in 0:size_ext-1
            if j<=shift
                a = j-shift
                b = source[1] + l1low*a
                if (b<0) b=0. end
                working_space[size_ext+j+1] = b - working_space[size_ext+j+1]
            elseif (j >= ssize + shift)
                a = j - (ssize-1+shift)
                b = source[ssize]
                if (b<0) b=0. end
                working_space[size_ext+j+1] = b - working_space[size_ext+j+1];
            else
                working_space[size_ext+j+1] = source[j-shift+1] - working_space[size_ext+j+1];
            end
        end
        for j in 0:size_ext-1
            if (working_space[size_ext+j+1] < 0) working_space[size_ext+j+1] = 0. end
        end
    end

    for i in 0:size_ext-1
        working_space[i + 6*size_ext+1] = working_space[i + size_ext + 1]
    end

    if markov
        for j in 0:size_ext-1
            working_space[2*size_ext+j+1] = working_space[size_ext + j + 1]
        end
        xmin::Int = 0
        xmax::Int = size_ext-1
        maxch::Float64 = 0.
        for i in 0:size_ext-1
            working_space[i+1] = 0.
            if (maxch < working_space[2*size_ext+i+1]) maxch = working_space[2*size_ext+i+1] end
            plocha += working_space[2*size_ext+i+1]
        end
        if (maxch == 0) error("'maxch' = 0") end

        nom::Float64 = 1.
        working_space[xmin+1]=1.
        for i in xmin:xmax-1
            nip::Float64 = working_space[2 * size_ext + i + 1] / maxch
            nim::Float64 = working_space[2 * size_ext + i + 1 + 1] / maxch
            sp::Float64 = 0.
            sm::Float64 = 0.
            for l in 1:averWindow
                a::Float64 = (i+l>xmax) ? (working_space[2 * size_ext + xmax + 1] / maxch) : (working_space[2 * size_ext + i + l + 1] / maxch)
                b::Float64 = a-nip
                a = (a+nip <= 0) ? 1 : sqrt(a+nip)

                b = exp(b/a)
                sp += b
                a = ((i-l+a)<xmin) ? (working_space[2 * size_ext + xmin + 1] / maxch) : (working_space[2 * size_ext + i - l + 1 + 1] / maxch)
                b = a-nim
                a = (a+nim<=0) ? 1 : sqrt(a+nim)
                b = exp(b/a)
                sm +=b
            end
            a = sp/sm
            a = working_space[i + 1 + 1] = working_space[i + 1] * a
            nom += a
        end
        for i in xmin:xmax
            working_space[i+1] = working_space[i+1] / nom
        end
        for j in 0:size_ext-1
            working_space[size_ext + j + 1] = working_space[j + 1] * plocha
        end
        for j in 0:size_ext-1
            working_space[2 * size_ext + j + 1] = working_space[size_ext + j + 1]
        end
        if backgroundRemove
            for i in 1:numberIterations
                for j in i:size_ext-i-1
                    a = working_space[size_ext + j + 1]
                    b = (working_space[size_ext + j - i + 1] + working_space[size_ext + j + i + 1]) / 2.0
                    if (b<a) a=b end
                    working_space[j+1] = a;
                end
                for j in i:size_ext-i-1
                    working_space[size_ext + j + 1] = working_space[j + 1]
                end
            end
            for j in 0:size_ext-1
                working_space[size_ext + j + 1] = working_space[2 * size_ext + j + 1] - working_space[size_ext + j + 1]
            end
         end
    end

    # deconvolution starts
    area::Float64 = 0
    lh_gold::Int = -1
    posit::Int = 0
    maximum::Float64 = 0
    maximum_decon::Float64 = 0
    lda::Float64 = 0
    ldb::Float64 = 0
    ldc::Float64 = 0
    # generate response vector
    for i in 0:size_ext-1
        lda = i - 3 * sigma
        lda = lda^2 / (2 * sigma^2)
        j = round(Int, 1000.0*exp(-lda), RoundDown)
        lda = j
        if (lda != 0) lh_gold = i+1. end
        working_space[i+1] = lda
        area += lda
        if lda > maximum
            maximum = lda
            posit = i
        end
    end
    # read source vector
    for i in 0:size_ext-1
        working_space[2*size_ext+i+1] = abs(working_space[size_ext + i + 1])
    end

    # create matrix at a*vector(transponse(b))
    i = lh_gold-1
    if (i > size_ext) i=size_ext end

    imin::Int = -i
    imax::Int = i
    for i in imin:imax
        lda = 0.
        jmin::Int = 0
        if (i<0) jmin = -i end
        jmax::Int = lh_gold-i
        if (jmax > lh_gold-1) jmax = lh_gold-1 end
        for j in jmin:jmax
            ldb = working_space[j+1]
            ldc = working_space[i+j+1]
            lda += ldb*ldc
        end
        working_space[size_ext + i - imin + 1] = lda;
    end

    # # create vector p
    i = lh_gold-1
    imin = -i
    imax = size_ext+i-1
    for i in imin:imax
        lda = 0.
        for j in 1:lh_gold-1
            ldb = working_space[j+1]
            k = i + j
            if ((k >= 0) && (k<size_ext))
                ldc = working_space[2*size_ext+k+1]
                lda += ldb*ldc
            end
        end
        working_space[4*size_ext+i-imin+1] = lda
    end

    # # move vector p
    for i in imin:imax
        working_space[2 * size_ext + i - imin + 1] = working_space[4 * size_ext + i - imin + 1]
    end
    # # initialization of resulting vector
    for i in 0:size_ext-1
        working_space[i+1] = 1.
    end
    # # START OF ITERATIONS
    for lindex = 0:deconIterations-1
        for i in 0:size_ext-1
            if ((abs(working_space[2*size_ext+i+1]) > 0.00001) && (abs(working_space[i+1]) > 0.00001))
                lda = 0.
                jmin = lh_gold-1
                if (jmin>i) jmin = i end
                jmin = -jmin
                jmax = lh_gold-1
                if (jmax>(size_ext-i-1)) jmax = size_ext-i-1 end
                for j in jmin:jmax
                    ldb = working_space[j+lh_gold+size_ext]
                    ldc = working_space[i+j+1]
                    lda += ldb*ldc
                end
                ldb = working_space[2*size_ext+i+1]
                if (lda!=0)
                    lda = ldb/lda
                else
                    lda = 0.
                end
                ldb = working_space[i+1]
                lda = lda*ldb
                working_space[3*size_ext+i+1] = lda
            end
        end
        for i in 0:size_ext-1
            working_space[i+1] = working_space[3*size_ext+i+1]
        end
    end

    # # shift resulting spectrum
    for i in 0:size_ext-1
        lda = working_space[i+1]
        j = i+posit
        j = j % size_ext
        working_space[size_ext+j+1] = lda
    end
    # # write back resulting spectrum
    maximum = 0.
    maximum_decon = 0.
    j = lh_gold-1
    for i in 0:size_ext-1
        if ((i >=shift) && (i < ssize+shift))
            working_space[i+1] = area*working_space[size_ext+i+j + 1]
            if (maximum_decon < working_space[i+1]) maximum_decon = working_space[i+1] end
            if (maximum < working_space[6 * size_ext + i + 1]) maximum = working_space[6 * size_ext + i + 1] end
        else
            working_space[i+1] = 0
        end
    end
    lda = 1.
    if (lda > threshold) lda = threshold end
    lda = lda/100.

    fMaxPeaks::Int = 10000
    fPosition = zeros(Float64, fMaxPeaks)
    fPositionX = zeros(Float64, fMaxPeaks)
    fPositionY = zeros(Float64, fMaxPeaks)
    fResolution = 1
    fHistogram = 0
    fNPeaks = 0
    # peak_index = 0# already defined at the start
    # # searching for peaks in deconvolved spectrum
    for i in 1:size_ext-1-1
        if ((working_space[i+1] > working_space[i]) && working_space[i+1] > working_space[i+1+1])
            if ((i >= shift) && i<ssize+shift)
                if ((working_space[i+1]>lda*maximum_decon) && (working_space[6*size_ext+i+1] > threshold*maximum / 100.))
                    a = 0.
                    b = 0.
                    for j in i-1:i+1
                        a += (j-shift)*working_space[j+1]
                        b += working_space[j+1]
                    end
                    a = a/b
                    if (a<0) a = 0. end
                    if (a>=ssize) a=ssize-1 end
                    if (peak_index == 0)
                        fPositionX[0+1] = a
                        peak_index = 1
                    else
                        priz = 0
                        j = 0
                        while ((priz == 0) && (j<peak_index))
                        # for j in 0:peak_index-1
                            if (working_space[6*size_ext+shift+round(Int, a, RoundDown)+1] > working_space[6*size_ext+shift+round(Int, fPositionX[j+1], RoundDown)+1])
                                priz = 1
                            end
                            j+=1
                        end
                        if (priz == 0)
                            if (j<fMaxPeaks) fPositionX[j+1] = a end
                        else
                            for k in peak_index:-1:j
                                if (k<fMaxPeaks) fPositionX[k+1] = fPositionX[k-1+1] end
                            end
                            fPositionX[j] = a
                        end
                        if (peak_index < fMaxPeaks) peak_index += 1 end
                    end
                end
            end
        end
    end

    # return working_space
    destVector = zeros(Float64, ssize)
    for i in 0:ssize-1
        destVector[i+1] = working_space[i + shift + 1]
    end
    fNPeaks = peak_index
    if (peak_index == fMaxPeaks) warn("Peak buffer full", prefix="Warning: SearchHighRes -> ") end

    return destVector, fPositionX[1:fNPeaks]
end

"""
    peakfinder(h::Histogram; <keyword arguments>)::Tuple{Histogram, Array{Float64, 1}} 

Returns a deconvoluted spectrum and an array of peak positions.
  
# Keywords
- `σ::Real=2.0`: The expected sigma of a peak in the spectrum. In units of bins. 
- `threshold::Real=10.0`: Threshold for being identified as a peak in the deconvoluted spectrum. A single bin is identified as a peak when its weight exceeds the `threshold` and the previous bin was not identified as an peak.
- `backgroundRemove::Bool=true`
- `deconIterations::Int=3`
- `markov::Bool=true`
- `averWindow::Int=3`

# Source
This function is basically a copy of `TSpectrum::SearchHighRes` from [ROOT](https://root.cern.ch/doc/master/classTSpectrum.html).

- M.A. Mariscotti: A method for identification of peaks in the presence of background and its application to spectrum analysis. NIM 50 (1967), 309-320.
- M. Morhac;, J. Kliman, V. Matouoek, M. Veselsky, I. Turzo.:Identification of peaks in multidimensional coincidence gamma-ray spectra. NIM, A443 (2000) 108-125.
- Z.K. Silagadze, A new algorithm for automatic photopeak searches. NIM A 376 (1996), 451.
"""
function peakfinder(h::Histogram; kw...)::Tuple{Histogram, Array{Float64, 1}} 
    destVector, fPositionX = peakfinder(h.weights; kw...)
    rh = Histogram(h.edges[1], destVector, :left)
    rh.weights = destVector

    mp = StatsBase.midpoints(h.edges[1]);
    fPositionX_idcs = Int.(round.(1 .+ fPositionX, RoundUp));
    weight_in_bin = (fPositionX .+ 2) .% fPositionX_idcs .- 0.5
    binvolumes = [StatsBase.binvolume(h, i) for i in fPositionX_idcs];

    fPositionX = mp[fPositionX_idcs] + binvolumes .* weight_in_bin - 0.5 * binvolumes

    return rh, fPositionX
end
