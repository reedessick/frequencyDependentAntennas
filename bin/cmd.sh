sanityCheck \
    -v \
    --Ndirc 5 \
    --nside 64 \
    -o test/ 

#---

plot_antennas \
    -V \
    --ex 1 0 0 \
    --ey 0 1 0 \
    --Ntrial 5 \
    --Nfreq 101 \
    --nside 128 \
    -o mollweide

plot_network_antennas \
    -v \
    -i H -i L \
    -o network/ -t HL \
    $(seq 0 1000 75000)

plot_network_antennas \
    -v \
    -i H -i L -i V \
    -o network/ -t HLV \
    $(seq 0 1000 75000)

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
    --theta-phi 90 0 --theta-phi 90 45 --theta-phi 90 90 --theta-phi 90 135 --theta-phi 90 180 --theta-phi 90 225 --theta-phi 90 270 --theta-phi 90 315 \
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
    --theta-phi 90 0 --theta-phi 90 45 --theta-phi 90 90 --theta-phi 90 135 --theta-phi 90 180 --theta-phi 90 225 --theta-phi 90 270 --theta-phi 90 315 \
    -o bode \
    --norm2zeroFreq -t normed 

#---

for phi in 0 7.5 15 22.5
do
    bode_antennas \
        -V \
        --ex 1 0 0 \
        --ey 0 1 0 \
        --Nfreq 1001 --fmin 0.01 \
        --min-mag 1e-2 \
        --theta-phi  0 $phi \
        --theta-phi 30 $phi \
        --theta-phi 60 $phi \
        --theta-phi 90 $phi \
        --series-legend \
        -o bode -t series-${phi}
done

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

for nside in 2 4 8 16
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
                --L ${L}e3 \
                -o snr -t ${L}km-nside${nside}-${N}x${N} 
        done
    done
done

#--- some basic localization

### values to inject
theta=40
phi=300
psi=0
iota=0
distance=50
tac=0

### values for sampler
Nprior=1000
Nsteps=1000
Nwalk=50
threads=10

### single IFO runs
for ifo in H L
do
    localize \
        --time -V \
        -o localize \
        -t $ifo \
        -i $ifo \
        --NpriorSteps $Nprior \
        --Nwalkers $Nwalk \
        --Nsteps $Nsteps \
        --threads $threads \
        --args-are-deg \
        $theta $phi $psi $iota $distance $tac

    plot_ensemble \
        -v 
        -o localize/ \
        -t $ifo \
        --theta $theta \
        --phi $phi \
        --psi $psi \
        --iota $iota \
        --distanceMpc $distance \
        --timeAtCoalescence $tac \
        --angles-are-deg \
        localize/localize_${ifo}.txt 
done

### 2-detector aLIGO network
localize \
    --time -V \
    -o localize \
    -t HL \
    -i H -i L \
    --NpriorSteps $Nprior \
    --Nwalkers $Nwalk \
    --Nsteps $Nsteps \
    --threads $threads \
    --args-are-deg \
    $theta $phi $psi $iota $distance $tac

plot_ensemble \
    -v 
    -o localize/ \
    -t HL \
    --theta $theta \
    --phi $phi \
    --psi $psi \
    --iota $iota \
    --distanceMpc $distance \
    --timeAtCoalescence $tac \
    --angles-are-deg \
    localize/localize_HL.txt 
