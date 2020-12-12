file=Monte_Carlo_ISING_1D.exe
if [ -f $file ]; then
	n=100
	tmax=2
	tmin=0.5
	path=../data/input.dat
	for i in $(seq $n); do
		echo -n 'T='
		echo "scale=3; $tmin+$i/$n*$tmax-$i/$n*$tmin" |bc
		echo "scale=3; $tmin+$i/$n*$tmax-$i/$n*$tmin" |bc > $path #temperature
		echo '50' >> $path #nspin
		echo '1.0' >> $path #J
		echo '0.02' >> $path #h
		echo '0' >> $path #metro
		echo '1' >> $path #blocks
		echo '1000' >> $path #steps per block (equilibration)
		echo '0' >> $path #restart
		echo '0' >> $path #instant
		echo '0' >> $path #vsTemperature
		./$file
		sed -i '8s/.*/1/' $path #attivo restart
		./$file
		if [ $i -eq 1 ] || [ $i -eq $n ]; then
			sed -i '9s/.*/1/' $path #scrivo l'ultima equilibrazione di T=0.515 e T=2
			./$file
		else
			./$file
		fi
		sed -i '6s/.*/10/' $path #veri blocchi e step della simulazione
		sed -i '7s/.*/10000/' $path
		sed -i '9s/.*/0/' $path
		sed -i '10s/.*/1/' $path
		./$file
	done
	
else
	echo 'Error: executable not found'
fi
