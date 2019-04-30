        comp := gfortran
        opt  := -c -Wall -Wtabs
        pattern:=*.f95
        source :=$(wildcard $(pattern))
        obj :=$(patsubst %.f95, %.o, $(source))

        main : $(obj)
	       $(comp) $^ -o $@

         %.o : %.f95
	       $(comp) $(opt) $< -o $@

       %.mod : %.f95
	       $(comp) $(opt) $<

    main.o : subprograms.mod

       result : main input
	        time ./$<  < input > output
       result-r : main input
		rm output
		make result
		cat output

        clean :
	 rm -f *.o *.mod main

        clean-all :
	 rm -f *.o *.mod main *.eps *.dat result

          git-s :
		git status
		git remote

        ifeq (git,$(firstword $(MAKECMDGOALS)))
        rep := $(wordlist 2,2,$(MAKECMDGOALS))
        m := $(wordlist 3,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
        $(eval $(a):;#)
        endif

          git :
		git add -A
		git commit -m "$(m)"
		git push -u $(rep) master

    git-clean :
		rm -rf .git

        ifeq (git-new,$(firstword $(MAKECMDGOALS)))
        new_rep := $(wordlist 2,2,$(MAKECMDGOALS))
        label := $(wordlist 3,3,$(MAKECMDGOALS))
        $(eval $(a):;#)
        endif

      git-new :
		make git-clean
		git init
		git remote add $(label) git@github.com:Paveloom/$(new_rep).git
		git add Makefile
		git commit -m "Стартовый make-файл."
		git push -u $(label) master
