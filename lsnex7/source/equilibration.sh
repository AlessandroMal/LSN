ex=equilibration.exe
inp=../data/configurations/input.dat
inpfile='../data/configurations/input.'$1
cp $inpfile $inp
sed -i "6s/.*/$2/" $inp #inserisci nblocks
sed -i "7s/.*/$3/" $inp #inserisci nsteps
echo 'Equilibrating system...'
./$ex
cp ../data/configurations/config.final ../data/configurations/config.0
for i in $(seq $4); do #inserisci quante volte voglio rilanciare
	./$ex
	cp ../data/configurations/config.final ../data/configurations/config.0
done
