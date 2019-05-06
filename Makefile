
      # Настройки компиляции программ
           comp := gfortran
            opt := -c -Wall -Wtabs
        pattern := f95
     allpattern := *.$(pattern)
     anypattern := %.$(pattern)
         source := $(wildcard $(allpattern))
            mod := $(patsubst $(anypattern), %.mod, $(source))
            obj := $(patsubst $(anypattern), %.o, $(source))

      # Блок правил для компиляции объектных файлов
      
        main : $(obj)
	       $(comp) $^ -o $@

         %.o : %.f95
	       $(comp) $(opt) $< -o $@

       %.mod : %.f95
	       $(comp) $(opt) $<

      # Блок правил-зависимостей (при необходимости)
        main.o : testm.mod

      # Блок правил для инициализации make-файла для сборки программы
      
           input :
	              touch input
      
          result : main input
	              time ./$<  < input > output
	        
        result-r : main input
		         rm output
		         make result
		         cat output

      # Блок правил для очистки директории
        clean :
	 rm -f $(obj) $(mod) main

        clean-all :
	 rm -f *.o *.mod main *.eps *.dat result



      # Блок правил для загрузки кода на Github

        # Правила для проверки статуса репозитория
          git-s :
		git status
		git remote

        # Правило для загрузки на Github с указанием метки репозитория и сообщения коммита (см. Readme)

        ifeq (git-r,$(firstword $(MAKECMDGOALS)))
        rep := $(wordlist 2,2,$(MAKECMDGOALS))
        m := $(wordlist 3,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
        $(eval $(rep):;#)
        $(eval $(m):;#)
        endif

          git-r :
		git add -A
		git commit -m "$(m)"
		git push -u $(rep) master
		
	  # Правило для загрузки на Github с указанием сообщения коммита, но без указания метки репозитория (см. Readme)
	  
        ifeq (git,$(firstword $(MAKECMDGOALS)))
        m := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
        $(eval $(m):;#)
        endif

          git :
		git add -A
		git commit -m "$(m)"
		git push -u origin master

        # Правило для удаления репозитория в текущей директории
          git-clean :
		rm -rf .git

       # Правило для подключения нового репозитория с указанием названия и метки
       # и загрузки в него стартового make-файла (см. Readme)

        ifeq (git-new-r,$(firstword $(MAKECMDGOALS)))
        new_rep := $(wordlist 2,2,$(MAKECMDGOALS))
        label := $(wordlist 3,3,$(MAKECMDGOALS))
        $(eval $(new_rep):;#)
        $(eval $(label):;# Успешно.)
        endif

         git-new-r :
			make git-clean
			git init
			git remote add $(label) git@github.com:Paveloom/$(new_rep).git
			git add Makefile
			git commit -m "Стартовый make-файл."
			git push -u $(label) master
			
       # Правило для подключения нового репозитория с указанием названия, но без указания метки,
       # и загрузки в него стартового make-файла (см. Readme)

        ifeq (git-new,$(firstword $(MAKECMDGOALS)))
        new_rep := $(wordlist 2,2,$(MAKECMDGOALS))
        $(eval $(new_rep):;# Успешно.)
        endif

         git-new :
			make git-clean
			git init
			git remote add origin git@github.com:Paveloom/$(new_rep).git
			git add Makefile
			git commit -m "Стартовый make-файл."
			git push -u origin master
			

