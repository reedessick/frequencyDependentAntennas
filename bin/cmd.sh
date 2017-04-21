plot_antennas \
    -V \
    --ex 1 0 0 \
    --ey 0 1 0 \
    --Ntrial 5 \
    --Nfreq 101 \
    --nside 128 \
    -o mollweide

#---

bode_antennas \
    -V \
    --ex 1 0 0 \
    --ey 0 1 0 \
    --Nfreq 1001 --fmin 0.01 \
    --min-mag 1e-4 \
    --theta-phi 0  0.0 --theta-phi 15  0.0 --theta-phi 30  0.0 --theta-phi 45  0.0 --theta-phi 60  0.0 --theta-phi 75  0.0 --theta-phi 90  0.0 \
    --theta-phi 0  7.5 --theta-phi 15  7.5 --theta-phi 30  7.5 --theta-phi 45  7.5 --theta-phi 60  7.5 --theta-phi 75  7.5 --theta-phi 90  7.5 \
    --theta-phi 0 15.0 --theta-phi 15 15.0 --theta-phi 30 15.0 --theta-phi 45 15.0 --theta-phi 60 15.0 --theta-phi 75 15.0 --theta-phi 90 15.0 \
    --theta-phi 0 22.5 --theta-phi 15 22.5 --theta-phi 30 22.5 --theta-phi 45 22.5 --theta-phi 60 22.5 --theta-phi 75 22.5 --theta-phi 90 22.5 \
    -o bode 

bode_antennas \
    -V \
    --ex 1 0 0 \
    --ey 0 1 0 \
    --Nfreq 1001 --fmin 0.01 \
    --min-mag 1e-4 \
    --theta-phi 0  0.0 --theta-phi 15  0.0 --theta-phi 30  0.0 --theta-phi 45  0.0 --theta-phi 60  0.0 --theta-phi 75  0.0 --theta-phi 90  0.0 \
    --theta-phi 0  7.5 --theta-phi 15  7.5 --theta-phi 30  7.5 --theta-phi 45  7.5 --theta-phi 60  7.5 --theta-phi 75  7.5 --theta-phi 90  7.5 \
    --theta-phi 0 15.0 --theta-phi 15 15.0 --theta-phi 30 15.0 --theta-phi 45 15.0 --theta-phi 60 15.0 --theta-phi 75 15.0 --theta-phi 90 15.0 \
    --theta-phi 0 22.5 --theta-phi 15 22.5 --theta-phi 30 22.5 --theta-phi 45 22.5 --theta-phi 60 22.5 --theta-phi 75 22.5 --theta-phi 90 22.5 \
    -o bode \
    --norm2zeroFreq -t normed 

#---

for L in 4 10 20 30 40 50
do

    plot_waveforms \
        -V \
        --ex 1 0 0 \
        --ey 0 1 0 \
        --theta-phi-psi 0 0 0 \
        --plot-tmin -0.25 \
        -o waveforms/ \
        -L ${L}e3 -t ${L}km
done

m1=10 
for m2 in 1 2 4 8 10 
do 
    for L in 4 10 20 30 40 50 
    do 
        plot_waveforms \
            -V \
            --ex 1 0 0 \
            --ey 0 10 0 \
            --theta-phi-psi 0 0 0 \
            --plot-tmin -0.25 \
            --m1 ${m1} \
            --m2 ${m2} \
            -L ${L}e3 \
            -o waveforms \
            -t ${m1}-${m2}_${L}km
    done 
done

#---

for L in 4 30 40 50 
do 
    snr_loss \
        -V \
        --theta-phi-psi-iota  0  0.0 0 0 \
        --theta-phi-psi-iota 30  0.0 0 0 \
        --theta-phi-psi-iota 60  0.0 0 0 \
        --theta-phi-psi-iota 90  0.0 0 0 \
        --theta-phi-psi-iota 30  7.5 0 0 \
        --theta-phi-psi-iota 60  7.5 0 0 \
        --theta-phi-psi-iota 90  7.5 0 0 \
        --theta-phi-psi-iota 30 15.0 0 0 \
        --theta-phi-psi-iota 60 15.0 0 0 \
        --theta-phi-psi-iota 90 15.  0 0 \
        --theta-phi-psi-iota 30 22.5 0 0 \
        --theta-phi-psi-iota 60 22.5 0 0 \
        --theta-phi-psi-iota 90 22.5 0 0 \
        -L ${L}e3 \
        -o snr -t ${L}km 
 done

#---

for N in 5 10 50 
do 
    for L in 4 30 40 50 
    do 
        snr_loss \
            -V \
            --theta-phi  0  0.0 \
            --theta-phi 30  0.0 \
            --theta-phi 60  0.0 \
            --theta-phi 90  0.0 \
            --theta-phi 30  7.5 \
            --theta-phi 60  7.5 \
            --theta-phi 90  7.5 \
            --theta-phi 30 15.0 \
            --theta-phi 60 15.0 \
            --theta-phi 90 15.0 \
            --theta-phi 30 22.5 \
            --theta-phi 60 22.5 \
            --theta-phi 90 22.5 \
            -L ${L}e3 \
            --N-psi $N \
            --N-iota $N \
            -o snr -t ${L}km-${N}x${N}
    done 
done

#---

for nside in 4 8 16
do
    for N in 5 10
    do
        for L in 4 30 40 50 
        do 
            snr_loss \
                -V \
                --heatmap \
                --nside $nside \
                --N-psi $N \
                --N-iota $N \
                --L $L \
                -o snr -t ${L}km-nside${nside}-${N}x${N} 
        done
    done
done
