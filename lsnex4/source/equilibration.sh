ex=equilibration.exe
inp=../data/configurations/input.dat
inpfile='../data/configurations/input.'$1
cp $inpfile $inp
rm ../data/configurations/old.0
sed -i "6s/.*/$2/" $inp #inserisci nblocks
sed -i "7s/.*/$3/" $inp #inserisci nsteps
sed -i "8s/.*/$3/" $inp #iprint
N=$( expr $4 + 1 )
echo 'Equilibrating system...'
for i in $(seq $N); do #inserisci quante volte voglio rilanciare
	./$ex
	cp ../data/configurations/config.final ../data/configurations/config.0
	cp ../data/configurations/old.final ../data/configurations/old.0
	mv ../data/measures/equil_temp.out ../data/measures/equil_temp$i.out
done
