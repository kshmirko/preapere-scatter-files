program main
	use mo_par_dls
	use mo_dls, only: key, keyel, ndp
	use arnt_hdrs
	implicit none
	integer, parameter	::	dp = kind(1.d0)
	
	integer(kind=kind(10))	::	I,J, ERR, unit, N, doy
	character(len=20)		::	site, date, time
	character(len=80)		::	errmsg
	character(len=15)		::	raw_row(COLUMNS_COUNT)
	real(kind=dp)			::	used_row(USED_COLUMNS_COUNT)
	logical 				::	error_read_line, &
								&error_parse_line
	N=10

	open(newunit=unit, file='data.csv', status='old', iostat=ERR,&
		& iomsg=errmsg)
	if (ERR.ne.0) then
    	
		print *, errmsg
    	stop
    endif

	! Skip header line
    read(unit, *)
	
	I=0
	error_parse_line = .FALSE.
	error_read_line = .FALSE.
mainloop: DO 
		! очищаем память под массив данных
		raw_row = ''
		! читаем строку из файла в массив строк
    	read(unit, *, iostat=ERR, iomsg=errmsg) &
			&(raw_row(J),J=1,COLUMNS_COUNT)
		
		! Если при чтении возникла ошибка, 
		if (ERR.ne.0) then
    		print '(A, I5)', 'Error=', ERR
			print '(A)', errmsg
			error_read_line = .TRUE.
			exit mainloop
		end if
		
		I=I+1 ! увеличиваем счетчик прочитанных строк
		
		! Для отладочных целей, в продакшене нужно убрать этот блок
		if (I .gt. 10) then
			exit mainloop
		end if
		
		!разбираем интересующие нас столбцы
		inner_loop: DO J=1, USED_COLUMNS_COUNT
			read(raw_row(indices(J)), *, iostat=ERR, iomsg=errmsg)&
			& used_row(J)
			
			
			if (ERR.ne.0) then
	    		print '(A, I5)', 'Error=', ERR
				print '(A)', errmsg
				error_parse_line=.TRUE.
				exit mainloop
			end if
			
		END DO inner_loop
		
		if (.not. error_parse_line) then
		end if
		
		!print '(4A10)', raw_row(indices(1):indices(4))
		!print '(4F7.2)', used_row(1:4)
    END DO mainloop

    close(unit)

	
	if (.not. (error_read_line .and. error_parse_line)) then
	    call DLS_read_input("input_sphrs.dat")
		call alloc_dls_array(key, keyel, 1)

		call optchar(ndp)

		call alloc_dls_array(key, keyel, 2)
	end if
    
	
end program main

