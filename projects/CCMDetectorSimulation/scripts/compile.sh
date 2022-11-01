
#/bin/bash

n=1

    for file in `find G4W*Sodium*.txt -mmin +2 -type f`
    do
	main=$file
	echo $main
	pat=G*_
	main2=${main##$pat}
#    echo $main2
	main3=tenkSodiumSimulation_${main2}
	echo $main3
	cat $file >> ${main3}
	rm $file
#    mv $main2 currentresults/
    done
    for file in `find G4W*Cobalt*.txt -mmin +2 -type f`
    do
	main=$file
	echo $main
	pat=G*_
	main2=${main##$pat}
#    echo $main2
	main3=tenkCobaltSimulation_${main2}
	echo $main3
	cat $file >> ${main3}
	rm $file
#    mv $main2 currentresults/
    done
    
    for file in `find G4W*Laser*.txt -mmin +2 -type f`
    do
	main=$file
	echo $main
	pat=G*_t
	main2=${main#$pat}
#    echo $main2
	main3=t${main2}
	echo $main3
	cat $file >> ${main3}
	rm $file
#    mv $main2 currentresults/
    done
