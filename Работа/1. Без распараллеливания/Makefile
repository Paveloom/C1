
     # Настройки компиляции программ
     
     comp := gfortran
     opt := -c -Wall -Wtabs
     pattern := f95
     allpattern := *.$(pattern)
     anypattern := %.$(pattern)
     source := $(wildcard $(allpattern))
     mod := $(patsubst $(anypattern), %.mod, $(source))
     obj := $(patsubst $(anypattern), %.o, $(source))
     
     # Заглушка на вывод сообщений указанными правилами
     # (без указания имён подавляет вывод со стороны make-файла у всех правил)
     .SILENT: 

     # Блок правил для компиляции объектных файлов
      
     main : $(obj)
	       $(comp) $^ -o $@

     %.o : $(anypattern)
	      $(comp) $(opt) $< -o $@

     %.mod : $(anypattern)
	        $(comp) $(opt) $<

     # Блок правил-зависимостей (при необходимости)
     main.o : subprograms.mod

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
			
			
			
     # Забирает последний коммит из master и копирует его в feature с тем же содержимым
     
     ifeq (git-ch,$(firstword $(MAKECMDGOALS)))
          line := $(shell git log master --oneline | head -n 1)
          hash := $(wordlist 1, 1, $(line))
          commit := $(wordlist 2, $(words $(line)), $(line))
          $(eval $(line):;@#)
          $(eval $(hash):;@#)
          $(eval $(commit):;@#)
     endif
     
     git-ch : 
	         git checkout features
	         echo $(hash)
	         echo $(commit)
	         git cherry-pick -n $(hash)
	         git add -A
	         git commit -m "$(commit)"
	         
	# Правило на работу в ветке feature
	
     git-f :
		 git add -A
		 git commit -e
		 git push --force-with-lease origin features

     git-fam :
	          git add -A
	          git commit --amend
	          git push --force-with-lease origin features  

