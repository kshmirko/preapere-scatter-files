module utils
	
contains
	subroutine FATAL(msg)
		implicit none
		character(*)	::	msg
		PRINT *, 100, msg
		
		100 FORMAT(T9, A)
		STOP
	end subroutine FATAL
end module	utils