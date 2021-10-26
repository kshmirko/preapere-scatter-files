program main
	use mo_par_dls
	use mo_dls, only: key, keyel, ndp
	use arnt_hdrs
	implicit none
	integer, parameter	::	dp = kind(1.d0)
	integer, parameter	::	MAX_RECS = 141
	integer(kind=kind(10))	::	I,J, ERR, unit, N, doy
	character(len=20)		::	site, date, time
	character(len=80)		::	errmsg
	real(kind=dp)			::	row(MAX_RECS)
	N=10

	
	
	open(newunit=unit, file='data.csv', status='old', iostat=ERR)
	if (ERR.ne.0) then
    	print *, 'Error opening file'
    	stop
    endif

	! Skip header line
    read(unit, *)
	
	I=0
	
mainloop: DO 
		dat = 0.0
    	read(unit, *, iostat=ERR, iomsg=errmsg) site, date, time, doy, &
		& (dat(J),j=1, MAX_RECS)
		if (ERR.ne.0) then
    		print '(A, I5)', 'Error=', ERR
			print '(A)', errmsg
			exit mainloop
		end if
		
		i=i+1
		if (I .gt. 10) then
			exit mainloop
		end if
		
		print '(3A15, 5F5.2)', site, date, time,&
		&row(AodExtTot_440nm:AodExtTot_1020nm)
    END DO mainloop

    close(unit)

    call DLS_read_input("input_sphrs.dat")
	call alloc_dls_array(key, keyel, 1)

	call optchar(ndp)

	call alloc_dls_array(key, keyel, 2)
	
end program main

