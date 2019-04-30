   comp:=gfortran
    opt:=-c -Wall -Wtabs
pattern:=*.f95
source :=$(wildcard $(pattern))
   obj :=$(patsubst %.f95, %.o, $(source))
      main : $(obj)
	     $(comp) $^ -o $@

       %.o : %.f95
		$(comp) $(opt) $< -o $@
     %.mod : %.f95
		$(comp) $(opt) $<

       result : main
	        time ./$< > output
       result-r : main
		rm output
		make result
		cat output
        clean :
	 rm -f *.o *.mod main

        clean-all :
	 rm -f *.o *.mod main *.eps *.dat result
	
