
     # Заглушка на вывод сообщений указанными правилами
     # (без указания имён подавляет вывод со стороны make-файла у всех правил)
     .SILENT: 

     # Настройки компиляции программ
     
     comp := mpiifort
     opt := -Ofast
     pattern := f90
     main_name := main
     
     # Директория с исходниками
     source_path := ./Исходники/
     
     # Правила-зависимости (при необходимости)
     $(source_path)$(main_name).o : $(source_path)subprograms.mod
     $(source_path)subprograms.mod : $(source_path)mpi_subprograms.mod
     
     # Распределение исходных файлов по переменным
     
     allpattern := $(source_path)*.$(pattern)
     anypattern := $(source_path)%.$(pattern)
     source := $(wildcard $(allpattern))
     mod := $(patsubst %.$(pattern), %.mod, $(filter-out %$(main_name).$(pattern), $(source)))
     obj_main := $(patsubst %.$(pattern), %.o, $(filter %$(main_name).$(pattern), $(source)))
     obj_mod := $(patsubst %.mod, %.o, $(mod))
     
     # Блок правил для компиляции объектных файлов
      
     $(main_name) : $(obj_main) $(obj_mod)
	               $(comp) $(opt) $^ -o $@

     $(source_path)%.mod : $(anypattern)
	                      $(comp) -c $(opt) -module $(source_path) $<
	                      if [ "./" != "$(source_path)" ]; then\
	                           mv -f $(patsubst $(source_path)%.mod, %.o, ./$@) $(source_path);\
	                      fi

     $(obj_main) : $(source_path)$(main_name).$(pattern)
	              $(comp) -c $(opt) $< -o $@
        
     # Блок правил для сборки программы
      
     input :
	        touch input
      
     result : $(main_name) input
	         time mpiexec -np 4 ./$< > output
	        
     result-r : input
		      rm output
		      make result
		      cat output

     # Блок правил для очистки директории
     
     clean :
	        rm -f $(obj_main) $(obj_mod) $(mod) $(main_name)

     clean-all :
	            rm -f *.o *.mod $(main_name) *.eps *.dat result



     # Блок правил для загрузки кода на Github

     # Правила для проверки статуса репозитория
          
     git-s :
		   git status
		   git remote

     # Правило для загрузки на Github с указанием метки репозитория и сообщения коммита (см. Readme)

     ifeq (git-r,$(firstword $(MAKECMDGOALS)))
          rep := $(wordlist 2,2,$(MAKECMDGOALS))
          $(eval $(rep):;@#)
     endif

     git-r :
		   git add -A
		   git commit -e
		   git push -u $(rep) master
		
	# Правило для загрузки на Github с указанием сообщения коммита, но без указания метки репозитория (см. Readme)

     git :
		 git add -A
		 git commit -e
		 git push -u origin master
		 
     # Правило для обновления последнего коммита до текущего состояния локального репозитория (см. Readme)
     # (без указания метки репозитория; использовать только при уверенности в безопасности)
		 
     git-am : 
	         git add -A
	         git commit --amend
	         git push --force-with-lease origin master
	         
     # Правило для обновления последнего коммита до текущего состояния локального репозитория (см. Readme)
     # (без указания метки репозитория; использовать только при уверенности в безопасности)

     ifeq (git-am-r,$(firstword $(MAKECMDGOALS)))
          label := $(wordlist 2,2,$(MAKECMDGOALS))
          $(eval $(label):;@#)
     endif
 
     git-am-r : 
	         git add -A
	         git commit --amend
	         git push --force-with-lease $(label) master

     # Правило для удаления репозитория в текущей директории
     git-clean :
		       rm -rf .git

     # Правило для подключения нового репозитория с указанием названия и метки
     # и загрузки в него стартового make-файла (см. Readme)

     ifeq (git-new-r,$(firstword $(MAKECMDGOALS)))
          new_rep := $(wordlist 2,2,$(MAKECMDGOALS))
          label := $(wordlist 3,3,$(MAKECMDGOALS))
          $(eval $(new_rep):;@#)
          $(eval $(label):;@#)
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
          $(eval $(new_rep):;@#)
     endif

     git-new :
			make git-clean
			git init
			git remote add origin git@github.com:Paveloom/$(new_rep).git
			git add Makefile
			git commit -m "Стартовый make-файл."
			git push -u origin master

