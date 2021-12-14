#====================================
#SOUND ANALYSIS OF ORGAN PIPES
#====================================
#Lin Yangchen
#December 2021
#for personal non-commercial use only



#====================================
#import wav files
#====================================

#if stereo, takes left channel
#if mono right channel, copies data to left channel and uses left channel

print(noquote(paste('importing wav files from', folder, 'folder')))

#numerical file names correspond to unique pipe IDs in pipedata
soundfiles <- list.files(pattern = '\\.wav$', ignore.case = T)
nclips <- length(soundfiles)

clipid <- str_remove(soundfiles, ".wav")
clipid <- str_remove(clipid, ".WAV")
clipid <- as.numeric(clipid)


sounddata    <- list(NULL)
maxtimestep  <- 0
minamplitude <- 0
maxamplitude <- 0
samprate     <- NULL
bitdepth     <- NULL
for(i in 1:nclips)
    {
        rawdata <- readWave(soundfiles[i])
        
        if(i == 1 &
        rawdata@samp.rate > downsamprate)
            {print(noquote(paste(
            'downsampling from', rawdata@samp.rate/1000,
            'kHz to', downsamprate/1000, 'kHz'
            )))}
        
        print(noquote(paste('pipe', clipid[i])))
        
        sampledata <- downsample(rawdata, downsamprate)
        
        if(length(sampledata@left) == 0)
            {
                sampledata@left <- sampledata@right
                sampledata@right <- numeric(0)
            }

        sounddata[[clipid[i]]] <- sampledata

        
        #no. of samples in longest recording
        maxtimestep <- max(maxtimestep, length(sounddata[[clipid[i]]]@left))

        #range of amplitude recorded across clips
        minamplitude <- min(minamplitude, min(sounddata[[clipid[i]]]@left))
        maxamplitude <- max(maxamplitude, max(sounddata[[clipid[i]]]@left))
        
        samprate <- c(samprate, sounddata[[clipid[i]]]@samp.rate)
        if(length(samprate) > 1 & var(samprate)!=0)
            {stop('clips have different sample rates')}

        bitdepth <- c(bitdepth, sounddata[[clipid[i]]]@bit)
        if(length(bitdepth) > 1 & var(bitdepth)!=0)
            {stop('clips have different bit depths')}

    }


maxtimestep <- maxtimestep - 1
samprate <- samprate[1]
bitdepth <- bitdepth[1]

save(sounddata, file = 'sounddata.RData')




#====================================
#waveforms
#====================================

print(noquote(paste('plotting pipe waveforms in', folder, 'folder')))

xmax <- maxtimestep/samprate
pdf(height = 5, width = 5, 'waveforms.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1))

    for(i in clipid)
        {
            print(noquote(paste('pipe', i)))
            
            x <- 0:(length(sounddata[[i]]@left) - 1)
            x <- x/samprate
            
            plot(x, sounddata[[i]]@left/maxamplitude
            , type = 'l', lwd = 0.25
            , xlim = c(0, xmax), xaxs = 'i'
            , ylim = c(minamplitude/maxamplitude, 1)
            , mgp = mgp, tck = tck, axes = F
            , xlab = 'time (s)', ylab = 'normalized sound pressure'
            )
            
            axis(1, mgp = mgp, tck = tck)
            axis(2, mgp = mgp, tck = tck)
            box(bty = 'L')
            
            text(xmax, 1,
            paste('pipe', i,
            samprate/1000, 'kHz', bitdepth, 'bit'),
            pos = 2)
        }



dev.off()
embedFonts('waveforms.pdf')





#====================================
#stationary spectra
#====================================

print(noquote(paste(
'RMS, fundamental frequencies and spectra/cepstra in',
folder, 'folder')))

samplerange <- timerange/xmax*maxtimestep

statspectra  <- list(NULL)
minspecamp   <- NULL
maxspecamp   <- NULL
fundamentals <- rep(NA, nrow(pipedata))
fundamps     <- rep(NA, nrow(pipedata))
rms          <- rep(NA, nrow(pipedata))

for(i in clipid)
    {
        print(noquote(paste('pipe', i, 'RMS and fundamental frequency')))
        
        subdata <- sounddata[[i]]@left[samplerange[1]:samplerange[2]]
        
        rms[i] <- sqrt(mean(subdata^2))
        
        hs <- spectrum(subdata
        #, span = #smoothing
        , plot = F)
        
        statspectra[[i]] <- hs
        
        minspecamp <- min(minspecamp, min(hs$spec))
        maxspecamp <- max(maxspecamp, max(hs$spec))
        
        #need to calculate only once
        if(i == clipid[1]){
        freqseq <- seq(0, samprate/2, length.out = length(hs$freq))}
        
        #determine fundamental frequency
        
        ref.freq <- pipedata$freq.eq.Hz[i]
        
        #look for peak
        freqinds <- which(freqseq > ref.freq*0.7 & freqseq < ref.freq*1.8)
        peak <- max(hs$spec[freqinds])
        fundamps[i] <- peak
        peakind <- which(hs$spec[freqinds] == peak)
        fundfreq <- freqseq[freqinds[peakind]]
        fundamentals[i] <- fundfreq
    }

save(statspectra, file = 'statspectra.RData')
save(fundamentals, file = 'fundamentals.RData')
save(fundamps, file = 'fundamps.RData')
save(rms, file = 'rms.RData')

#plot individual stationary spectra
pdf(height = 5, width = 5, 'stationary_spectra.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1)
)

statcepstra <- list(NULL)
mincepsamp   <- NULL
maxcepsamp   <- NULL
for(i in clipid)
    {
        print(noquote(paste('pipe', i, 'spectrum')))
        
        #plot frequency cut-off at specified harmonic
        if(!is.na(hcutoff))
        {
            freqcutoff <- fundamentals[i]*hcutoff
        } else
        {
            freqcutoff <- samprate/2
        }
        freqcutind <- which(freqseq <= freqcutoff)
        
        ylim = c(log10(minspecamp), log10(maxspecamp))*10
        
        plot(freqseq[freqcutind], log10(statspectra[[i]]$spec)[freqcutind]*10
        , type = 'l', lwd = 0.5
        , xlim = c(0, freqcutoff)
        , ylim = ylim
        , xaxs = 'i'
        , mgp = mgp, tck = tck, axes = F
        , xlab = 'frequency (Hz)', ylab = 'dB'
        )
        
        axis(1, mgp = mgp, tck = tck)
        axis(2, mgp = mgp, tck = tck)
        box(bty = 'L')
        
        text(fundamentals[i],
        log10(statspectra[[i]]$spec)[which(freqseq == fundamentals[i])]*10,
        sprintf('%0.2f Hz', round(fundamentals[i], digits = 2)),
        #labels = paste(round(fundamentals[i], digits = 2), 'Hz'),
        pos = 4
        )
        
        text(freqcutoff, max(ylim),
        paste('pipe', i, 'to harmonic', hcutoff - 1)
        , pos = 2)
        
        
        
        #cepstrum
        
        cepstrum <- spectrum(log10(statspectra[[i]]$spec)*10
        #, span = #smoothing
        , plot = F
        )
        
        statcepstra[[i]] <- cepstrum
        
        mincepsamp <- min(mincepsamp, min(cepstrum$spec))
        maxcepsamp <- max(maxcepsamp, max(cepstrum$spec))

        #need to calculate only once
        if(i == clipid[1]){
        maxquef <- length(freqseq)/max(freqseq)/2
        quefseq <- seq(0, maxquef, length.out = length(cepstrum$freq))
        }
        
        
    
    }

dev.off()
embedFonts('stationary_spectra.pdf')

save(statcepstra, file = 'statcepstra.RData')




function(){

#plot cepstra

pdf(height = 5, width = 5, 'stationary_cepstra.pdf', family = font)
par(pty = 's', mar = c(3.1,3.1,1,1))

for(i in clipid)
    {
        print(noquote(paste('pipe', i, 'cepstrum')))

        #ylim = c(mincepsamp, maxcepsamp)
        ylim <- c(0, 1000) #manual adjustment
        
        plot(quefseq, statcepstra[[i]]$spec
        , type = 'l'
        , xaxs = 'i'
        , xlim = c(0.02, 0.5), ylim = ylim
        , mgp = mgp, tck = tck, axes = F
        , xlab = 'quefrency (s)', ylab = expression(dB^2)
        )

        axis(1, mgp = mgp, tck = tck)
        axis(2, mgp = mgp, tck = tck)
        box(bty = 'L')
        
        text(max(quefseq), max(ylim),
        paste('pipe', i)
        , pos = 2)

    }

dev.off()
embedFonts('stationary_cepstra.pdf')

} #end of comment-out









#====================================
#flutter detection
#====================================

print(noquote(paste('detecting pipe flutter in',
folder, 'folder')))

pdf(height = 5, width = 5, 'flutter.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1))

for(i in clipid)
    {
        print(noquote(paste('pipe', i)))
        
        #freq cut-off below/above fundamental
        bandleft  <- fundamentals[i] - flutterband
        bandright <- fundamentals[i] + flutterband
        freqcutind <- which(freqseq > bandleft &
        freqseq < bandright)
        
        yvals <- log10(statspectra[[i]]$spec)[freqcutind]*10
        plot(freqseq[freqcutind], yvals
        , type = 'l', lwd = 0.5
        , xaxs = 'i'
        , mgp = mgp, tck = tck, axes = F
        , xlab = 'frequency (Hz)', ylab = 'dB'
        )
        
        axis(1, mgp = mgp, tck = tck)
        axis(2, mgp = mgp, tck = tck)
        box(bty = 'L')
        
        text(max(freqseq[freqcutind]), max(yvals),
        paste('pipe', i)
        , pos = 2)
    
    }
    
dev.off()
embedFonts('flutter.pdf')








#====================================
#voicing
#====================================


#3D stack of stationary power spectra of pipes in each stop
print(noquote(paste('plotting stacked voicing spectra/cepstra in',
folder, 'folder')))

pdf(height = 5, width = 5, 'voicing_spectra.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1))

stops <- unique(pipedata$stop[clipid])
spec3d <- list(NULL)
ceps3d <- list(NULL)
for(i in stops)
    {
        print(noquote(paste(i, 'spectra')))
        
        pipeid <- clipid[pipedata$stop[clipid] %in% i]
        
        #use only the specified no. of notes from the bottom of the keyboard
        keys <- pipedata$key[pipeid]
        keys <- keys[which(keys <= keycut)]

        #frequency cut-off at specified harmonic
        if(is.na(stackfreqcut))
            {
                stackfreqcut <- samprate/2
            }
        stackfreqcutind <- which(freqseq <= stackfreqcut)

        numfreqbin <- length(freqseq[stackfreqcutind])
        numquefbin <- length(quefseq)
        
        spec3d[[i]] <- matrix(NA, ncol = max(keys), nrow = numfreqbin)
        ceps3d[[i]] <- matrix(NA, ncol = max(keys), nrow = numquefbin)
        
        for(j in 1:length(keys))
            {
                #retrieve spectrum
                spec3d[[i]][,keys[j]] <- statspectra[[pipeid[j]]]$spec[stackfreqcutind]
                
                #retrieve cepstrum
                ceps3d[[i]][,keys[j]] <- statcepstra[[pipeid[j]]]$spec    
            }
        
        #plot spectra
        
        image(freqseq[stackfreqcutind], 1:max(keys)
        , log10(spec3d[[i]])*10
        , col = timPalette(10000)
        , mgp = mgp, tck = -0.01, axes = F
        , xlab = 'frequency (Hz)', ylab = 'note'
        , main = i
        )
        
        axis(1, mgp = mgp, tck = -0.01)
        axis(2, mgp = mgp, tck = -0.01)
        box(bty = 'L')

    }


dev.off()
embedFonts('voicing_spectra.pdf')



function(){

#plot cepstra

pdf(height = 5, width = 5, 'voicing_cepstra.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1))

for(i in stops)
    {

        print(noquote(paste(i, 'cepstra')))

        pipeid <- clipid[pipedata$stop[clipid] %in% i]
        
        image(quefseq, 1:max(keys)
        , log10(ceps3d[[i]])
        , col = timPalette(10000)
        , mgp = mgp, tck = -0.01, axes = F
        , xlab = 'quefrency (s)', ylab = 'pipe no.'
        , main = i
        )

        axis(1, mgp = mgp, tck = -0.01)
        axis(2, mgp = mgp, tck = -0.01)
        box(bty = 'L')

    }

dev.off()
embedFonts('voicing_cepstra.pdf')

} #end of comment-out







#====================================
#tuning and volume
#====================================

print(noquote(paste('analyzing temperament of each stop in',
folder, 'folder')))

#use only the specified no. of notes from the bottom of the keyboard
exclind <- NULL
for(i in stops)
    {
        stopind <- which(pipedata$stop == i)
        exclind <- c(exclind, stopind[1:keycut])
    }
truncated <- fundamentals
truncated[-exclind] <- NA
pipedata$actual.freq <- truncated

#semitone difference from reference frequencies
pipedata$out.of.tune <- 1200*log2(pipedata$actual.freq/pipedata$freq.eq.Hz)/100

#normalize to A = 440 Hz
pipedata$out.of.tune <- pipedata$out.of.tune - mean(pipedata$out.of.tune, na.rm = T)


#compare with equal temperament

#root mean square error
tuningrmse <- sqrt(mean(pipedata$out.of.tune^2, na.rm = T))


pipedata$rms <- rms/max(rms, na.rm = T)
pipedata$fundamp.norm <- fundamps/max(fundamps, na.rm = T)


pdf(width = 9, height = 5, 'temperament_vol.pdf', family = font)

function(){

print(
xyplot(freq.eq.Hz ~ key | stop, data = pipedata
, type = 'l', col = 'black'
, aspect = 1
, xlab = 'note', ylab = 'frequency (Hz)'
, par.settings = list(strip.background = list(col = 'transparent'))
) +
xyplot(actual.freq ~ key | stop, data = pipedata
, col = 'red'
)
)

} #end of comment-out


#how many semitones out of tune
print(
barchart(out.of.tune ~ key | stop, data = pipedata
, col = 'black'
, aspect = 1
, xlab = 'note', ylab = 'tuning error (semitone)'
, scales = list(x = list(draw = F))
, par.settings = list(strip.background = list(col = 'transparent'))
, origin = 0, horizontal = F, reference = F
)
)


#tuning error distribution
print(
densityplot(~ out.of.tune | stop, data = pipedata
, type = 'density'
, col = 'black'
, aspect = 1
, xlab = 'tuning error (semitone)', ylab = 'probability density'
, par.settings = list(strip.background = list(col = 'transparent'))
)
)



#stationary RMS volume

print(
barchart(rms ~ key | stop, data = pipedata
, col = 'black'
, aspect = 1
, xlab = 'note', ylab = 'normalized RMS sound pressure'
, scales = list(x = list(draw = F))
, par.settings = list(strip.background = list(col = 'transparent'))
, origin = 0, horizontal = F, reference = F
)
)



#stationary fundamental volume

print(
barchart(fundamp.norm ~ key | stop, data = pipedata
, col = 'black'
, aspect = 1
, xlab = 'note', ylab = 'normalized fundamental sound pressure'
, scales = list(x = list(draw = F))
, par.settings = list(strip.background = list(col = 'transparent'))
, origin = 0, horizontal = F, reference = F
)
)

dev.off()
embedFonts('temperament_vol.pdf')








#====================================
#attack transients
#====================================

print(noquote(paste('plotting attack transient spectrograms in',
folder, 'folder')))

ind1 <- samprate*attacktime[1]
ind2 <- samprate*attacktime[2]


pdf(height = 5, width = 5, 'attack_spectrograms.pdf', family = font)
par(pty = 's', mar = c(3,3,1,1))

#power spectrogram
spectrograms <- list(NULL)
for(i in clipid)
    {
        print(noquote(paste('pipe', i)))
        
        spectrogram <- powspec(sounddata[[i]]@left[ind1:ind2], sr = samprate)
        spectrograms[[i]] <- spectrogram
        
        #frequency cutoff at specified harmonic
        afreqcutoff <- fundamentals[i]*ahcutoff
        cutoffrow <- min(1, afreqcutoff/(samprate/2))
        cutoffrow <- floor(nrow(spectrogram)*cutoffrow)
        
        image(seq(attacktime[1], attacktime[2], length.out = ncol(spectrogram))
        , seq(0, samprate/2, length.out = nrow(spectrogram))[1:cutoffrow]
        , log10(t(spectrogram))[,1:cutoffrow]
        , col = timPalette(10000)
        , mgp = mgp, tck = -0.01, axes = F
        , xlab = 'time (s)', ylab = 'frequency (Hz)'
        )

        axis(1, mgp = mgp, tck = -0.01)
        axis(2, mgp = mgp, tck = -0.01)
        box(bty = 'L')
        
        mtext(paste('pipe', i, 'to harmonic', ahcutoff - 1))
        
    }


dev.off()
embedFonts('attack_spectrograms.pdf')


save(spectrograms, file = 'spectrograms.RData')










save(pipedata, file = 'pipedata.RData')