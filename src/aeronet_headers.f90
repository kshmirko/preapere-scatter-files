module arnt_hdrs
! Наименование столбцов данныж
	integer, parameter	::	AodExtTot_440nm = 1
	integer, parameter	::	AodExtTot_675nm = 2
	integer, parameter	::	AodExtTot_870nm = 3
	integer, parameter	::	AodExtTot_1020nm =4
	integer, parameter	::	AngExpTot		= 5
	integer, parameter	::	AodAbsTot_440nm = 6
	integer, parameter	::	AodAbsTot_675nm = 7
	integer, parameter	::	AodAbsTot_870nm = 8
	integer, parameter	::	AodAbsTot_1020nm =9
	integer, parameter	::	RealMidx_440nm	= 10
	integer, parameter	::	RealMidx_675nm	= 11
	integer, parameter	::	RealMidx_870nm	= 12
	integer, parameter	::	RealMidx_1020nm	= 13
	integer, parameter	::	ImagMidx_440nm	= 14
	integer, parameter	::	ImagMidx_675nm	= 15
	integer, parameter	::	ImagMidx_870nm	= 16
	integer, parameter	::	ImagMidx_1020nm	= 17
	integer, parameter	::	Sphericity		= 18
	integer, parameter	::	R01				= 19
	integer, parameter	::	R02				= 20
	integer, parameter	::	R03				= 21
	integer, parameter	::	R04				= 22
	integer, parameter	::	R05				= 23
	integer, parameter	::	R06				= 24
	integer, parameter	::	R07				= 25
	integer, parameter	::	R08				= 26
	integer, parameter	::	R09				= 27
	integer, parameter	::	R10				= 28
	integer, parameter	::	R11				= 29
	integer, parameter	::	R12				= 30
	integer, parameter	::	R13				= 31
	integer, parameter	::	R14				= 32
	integer, parameter	::	R15				= 33
	integer, parameter	::	R16				= 34
	integer, parameter	::	R17				= 35
	integer, parameter	::	R18				= 36
	integer, parameter	::	R19				= 37
	integer, parameter	::	R20				= 38
	integer, parameter	::	R21				= 39
	integer, parameter	::	R22				= 40
	integer, parameter	::	VolCT			= 41
	integer, parameter	::	REffT			= 42
	integer, parameter	::	VMRT			= 43
	integer, parameter	::	StdT			= 44
	integer, parameter	::	VolCF			= 45
	integer, parameter	::	REffF			= 46
	integer, parameter	::	VMRF			= 47
	integer, parameter	::	StdF			= 48
	integer, parameter	::	VolCC			= 49
	integer, parameter	::	REffC			= 50
	integer, parameter	::	VMRC			= 51
	integer, parameter	::	StdC			= 52
	integer, parameter	::	LR_440nm		= 53
	integer, parameter	::	LR_675nm		= 54
	integer, parameter	::	LR_870nm		= 55
	integer, parameter	::	LR_1020nm		= 56
	integer, parameter	::	DLR_440nm		= 57
	integer, parameter	::	DLR_675nm		= 58
	integer, parameter	::	DLR_870nm		= 59
	integer, parameter	::	DLR_1020nm		= 60
	integer, parameter	::	SufAlbedo_440nm = 61
	integer, parameter	::	SufAlbedo_675nm = 62
	integer, parameter	::	SufAlbedo_870nm = 63
	integer, parameter	::	SufAlbedo_1020nm = 64
	integer, parameter	::	USED_COLUMNS_COUNT = 64
	integer, parameter	::	COLUMNS_COUNT 	= 155
    ! Номера столбцов данных из файла AERONET
	integer, parameter	::	indices(USED_COLUMNS_COUNT)		= &
	&[ 11, 12, 13, 14, 23, 28, 29, 30, 31, 33, 34, 35, 36, 37,&
	& 38, 39, 40, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, &
	& 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, &
	& 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 113, 114, 115,  &
	& 116, 117, 118, 119, 120, 142, 143, 144, 145 ]
	
end module arnt_hdrs
