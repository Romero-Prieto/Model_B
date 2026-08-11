tempfile   lISt
tempfile   temp
tempfile   temp2
tempfile   tempBR
tempfile   tempIR

local      pATh         = "/Users/lshjr3/Documents/DHS"
local      DHS          = "AF71 AL51 AL71 AM42 AM54 AM61 AM72 AO71 AO81 AZ52 BD31 BD3A BD41 BD4J BD51 BD61 BD72 BD7R BD81 BF21 BF31 BF43 BF62 BF81 BJ31 BJ41 BJ51 BJ61 BJ71 BO01 BO31 BO3B BO41 BO51 BR01 BR21 BR31 BU01 BU61 BU71 CD51 CD61 CD81 CF31 CG51 CG61 CI35 CI3A CI62 CI81 CM22 CM31 CM44 CM61 CM71 CO01 CO22 CO31 CO41 CO53 CO61 CO72 DR01 DR21 DR32 DR41 DR4B DR52 DR5A DR61 DR6A EC01 EG01 EG21 EG33 EG42 EG4A EG51 EG5A EG61 ES01 ET41 ET51 ET61 ET71 ET81 ET8A GA41 GA61 GA71 GH02 GH31 GH41 GH4B GH5A GH72 GH8C GM61 GM81 GN41 GN52 GN62 GN71 GU01 GU34 GU41 GU71 GY5I HN52 HN62 HT31 HT42 HT52 HT61 HT71 IA23 IA42 IA52 IA74 IA7E IAAP IAAR IAAS IABH IADL IAGJ IAGO IAHP IAHR IAJM IAKA IAKE IAMG IAMH IAMN IAMP IAMZ IANA IAOR IAPJ IARJ IASK IATN IATR IAUP IAWB ID01 ID21 ID31 ID3A ID42 ID51 ID63 ID71 JO21 JO31 JO42 JO51 JO61 JO6C JO74 JO81 KE03 KE33 KE3A KE42 KE52 KE72 KE8C KH42 KH51 KH61 KH73 KH82 KK31 KK42 KM32 KM61 KY31 KY61 LB01 LB51 LB6A LB7A LK02 LS41 LS61 LS71 LS81 MA01 MA21 MA43 MB53 MD21 MD31 MD42 MD51 MD81 ML01 ML32 ML41 ML53 ML6A ML7A ML8A MM71 MR71 MV52 MV71 MW22 MW41 MW4E MW61 MW7A MW81 MX01 MZ31 MZ41 MZ62 MZ81 NC31 NC41 NG21 NG4B NG53 NG6A NG7B NG8B NI22 NI31 NI51 NI61 NM21 NM41 NM51 NM61 NP31 NP41 NP51 NP61 NP7H NP82 OS01 PE01 PE21 PE31 PE41 PE51 PE5I PE61 PE6A PE6I PG71 PH31 PH3B PH41 PH52 PH61 PH71 PH82 PK21 PK53 PK61 PK71 PY21 RW21 RW41 RW53 RW5A RW61 RW70 RW81 SD02 SL51 SL61 SL7A SN02 SN21 SN32 SN4A SN61 SN6D SN6R SN71 SN7A SN7I SN7Z SN81 SN8B SN8S SNG0 ST51 SZ51 TD31 TD41 TD71 TG01 TG31 TG61 TH01 TJ61 TJ72 TJ81 TL61 TL71 TN02 TR31 TR41 TR4A TR51 TR62 TR71 TT01 TZ21 TZ3A TZ41 TZ4I TZ63 TZ7B TZ82 UA51 UG01 UG33 UG41 UG52 UG61 UG7B UZ31 VN31 VN41 YE21 YE61 ZA31 ZA71 ZM21 ZM31 ZM42 ZM51 ZM61 ZM71 ZM81 ZW01 ZW31 ZW42 ZW52 ZW62 ZW72"

clear
generate   survey       = ""
save      `lISt', replace

local      vARsBR       = "caseid v000 v001 v002 v003 v005 v008 v016 v011 v012 v022 v023 v024 v025 bidx b3 b4 b6 b17 v101 v102" /*Variables from Birth Recode.*/
local      vARsIR       = "caseid v000 v001 v002 v003 v005 v008 v016 v011 v012 v022 v023 v024 v025 v018 v101 v102 vcal_1" /*Variables from Individual Recode.*/
use       `vARsBR' in 1/1 using "`pATh'/Survey/BDBR81FL.DTA", clear      /*Most updated DHS. This file is created to make room for all variales.*/
drop       in 1                                                          /*Makes an empty file, just with variable names.*/
save      `tempBR', replace                                              /*Saves the file as a temporary file.*/

use       `vARsIR' in 1/1 using "`pATh'/Survey/BDIR81FL.DTA", clear      /*Follows the same approach with the IR file.*/
drop       in 1
save      `tempIR', replace


foreach svy of local DHS {
	cls
	dis      "`svy'"
	
	local      name           = substr("`svy'",1,2) + "BR" + substr("`svy'",3,4) + "FL.DTA"
	use        caseid v0* b* using "`pATh'/Survey/`name'", clear         /*Imports all variables starting by v0 and b.*/
	duplicates drop                                                      /*DHS samples do not have duplicates. This is only to show the appropriate place for this cleaning step, e.g. DR21 has 10 duplicates.*/
	merge 1:1  caseid v002 v012 bidx using `tempBR', nogenerate noreport nolabel keep(master) /*Makes room for the variables that are in the list but not in this specific sample.*/
	keep      `vARsBR'                                                   /*Keeps just the variables of the list.*/
	generate   mother         = 1                                        /*Creates an identifier of women in the BR. They should be mothers or mothers to be.*/
	save      `temp', replace                                            /*Saves date as a temporary file to be used later.*/

	local      name           = substr("`svy'",1,2) + "IR" + substr("`svy'",3,4) + "FL.DTA"
	use        caseid v0* b* using "`pATh'/Survey/`name'", clear         /*Imports all variables starting by v and b. IR includes all women aged 15–49, regardless of maternity.*/
	duplicates drop                                                      /*DHS samples do not have duplicates. This is only to show the appropriate place for this cleaning step.*/
	merge 1:1  caseid v002 v012 using `tempIR', nogenerate noreport nolabel keep(master) /*Brings the name of those variables not in the sample.*/
	keep      `vARsIR'                                                   /*Keeps just the variables of the list.*/
	merge 1:m  caseid v001 v002 v003 v008 v016 v008 v012 v023 v024 v025 using `temp', nogenerate /*Combines the IR and the BR into a single file.*/
	
	generate   fILe           = "`svy'"                                  /*Generates a indentifier for each survey.*/
	sort       caseid bidx                                               /*Assumes a natural sort. The code of the mother and the code of each birth.*/
	bysort     caseid: generate   k = _n                                 /*Creates an index. One number per record of a woman.*/
	bysort     caseid: generate   K = _N                                 /*Creates an index. Total number of records per woman.*/
	generate   cluster        = v001                                     /*DHS ids for cluster and stratum.*/
	
	sum        v023 
	replace    v023           = v022                                            if r(sd) == 0
	generate   stratum        = v023
		
	generate   w              = 1                                               if k     == 1
	recode     w           (. = 0)	
	sort       stratum cluster caseid k 
	generate   woman          = sum(w)
	drop       w
	
	capture    decode     v024, generate(state)
	capture    generate   state = ""
	replace    state          = "SS" + state + "SS"                             if state != ""
	
	generate   DHS            = substr(v000,1,2)
	merge m:1  DHS using "`pATh'/_codesDHS.dta", nogenerate noreport keep(master match) keepusing(country)
	drop       DHS

	generate   DOB            = mdy(v011 - floor((v011 - 1)/12)*12,1,floor((v011 - 1)/12) + 1900) /*Generates the date of birth of the woman/mother.*/
	recode     v016        (. = 1)                                       /*Exact date only available for the most recent surveys. If not avialable, the first of the month is assumed.*/
	generate   interview      = mdy(v008 - floor((v008 - 1)/12)*12,v016,floor((v008 - 1)/12) + 1900) /*Generates the date of the interview.*/
	replace    interview      = mdy(v008 + 1 - floor((v008 - 1)/12)*12,1,floor((v008 - 1)/12) + 1900) - 1        if interview == . /*Fix the date of the interview, for months with less days than those reported.*/
	
	/*Date of Birth - Birth Histories*/
	sort       b3 b17
	if b17[1] == . {
		/*The exact date is only avaliable for the most recent surveys. If not reported, a random day of the month is assumed. Identifies the limits of this random date.*/
		generate   B_min          = mdy(b3 - floor((b3 - 1)/12)*12,1,floor((b3 - 1)/12) + 1900)          if b3    != .            /*The first day of the reported month is the lower limit.*/
		generate   B_max          = mdy(b3 + 1 - floor(b3/12)*12,1,floor(b3/12) + 1900)                  if b3    != .            /*The first day of the following month is the upper limit.*/
		}
	else {
		generate   B_min          = mdy(b3 - floor((b3 - 1)/12)*12,b17,floor((b3 - 1)/12) + 1900)        if b3    != .            /*Calculates the date of birth.*/
		replace    B_min          = mdy(b3 + 1 - floor((b3 - 1)/12)*12,1,floor((b3 - 1)/12) + 1900) - 1  if B_min == . & b3  != . /*Calculates again, if incorrect day of the month.*/ 
		generate   B_max          = B_min                                                                                         /*Because exact dates are available.*/
		}
	replace    B_max          = max(min(B_max,interview),B_min)	                         if B_min != .                            /*Adjusts B_max postdating the day of interview.*/
		
	format     %tdDD/NN/CCYY interview B_* DOB		                     /*Gives date format to the date variables.*/
	rename     b4 sex                                                    /*Identifies the sex and age of the child.*/
	rename     v012 age
	generate   W              = v005/1000000                             /*Identifies the sampling weights, rounded and *10^6. Useful when decimals are not available, not the case.*/
	
	/*Ages at death (in days) - Birth Histories*/
	generate   D_min          = b6 - 100                                                 if b6    != .   & b6 <= 200 /*The report is in days within the first month of life.*/
	generate   D_max          = D_min + 1                                                if b6    != .   & b6 <= 200 /*Assumes a plausible maximum of one additional day.*/
	replace    D_min          = 0                                                        if b6    == 198 | b6 == 199 /*If days were reported but a number was not provided (rare).*/ 
	replace    D_max          = 365.25/12                                                if b6    == 198 | b6 == 199 /*Max and min bound the first month of life.*/	
	replace    D_min          = (b6 - 200)*365.25/12                                     if b6    != .   & b6 >= 200 & b6  < 300 /*The report is 2-24 months.*/
	replace    D_max          = (b6 - 200 + 1)*365.25/12                                 if b6    != .   & b6 >= 200 & b6  < 300 /*Assumes a maximum of one additional month.*/
	replace    D_min          = 0                                                        if b6    == 298 | b6 == 299 /*If months were reported but a number was not provided (rare).*/
	replace    D_max          = 24*365.25/12                                             if b6    == 298 | b6 == 299 /*Max and min bound the first 2 years of life.*/
	replace    D_min          = (b6 - 300)*365.25                                        if b6    != .   & b6 >= 300 /*The report is in years after the second birthday.*/
	replace    D_max          = (b6 - 300 + 1)*365.25                                    if b6    != .   & b6 >= 300 /*Assumes a plausible maximum of one additional year.*/
	replace    D_min          = 0                                                        if b6    == 398 | b6 == 399 /*If years were reported but a number was not provided (rare).*/
	replace    D_max          = max(year(interview) - year(B_min),0)*365.25              if b6    == 398 | b6 == 399 /*Age at death could be from 0 to the age at interview.*/
	replace    D_max          = max(max(min(B_max + D_max,interview) - B_max,0),D_min)   if D_min != .               /*Adjusts max ages at death postdating the day of interview.*/
	
	sort       caseid bidx                                                                                  /*Assumes a natural sort.*/
	generate   temp           = 1                                               if B_min != .               /*Creates a temp variable to indicate a child is born.*/
	bysort     caseid: egen       Born  = sum(temp)                                                         /*Creates a Summary of children ever born.*/ 
	replace    Born           = .                                               if k     != 1               /*Constraints the information to be one total per woman.*/
	recode     temp        (1 = .)                                              if b6    != .               /*Constratints the temp variable to indicate a child is alive.*/
	bysort     caseid: egen       Alive = sum(temp)                                                         /*Creates a Summary of children alive.*/ 
	replace    Alive          = .                                               if k     != 1               /*Constraints the information to be one total per woman.*/

	rename     v000 survey
	rename     v024 Region
	rename     v025 Urban
	
	keep       fILe survey country Region Urban state stratum cluster woman mother k K age sex interview B_* D_* W bidx caseid DOB Alive Born
	order      fILe survey country Region Urban state stratum cluster woman mother k K age sex interview B_* D_* W bidx caseid DOB Alive Born
	sort       fILe stratum cluster woman k
	label drop _all

	export     delimited using "`pATh'/OuTPuT/`svy'.csv", replace
	save      `temp2', replace

	if "`svy'" == "IA7E" | "`svy'" == "IA52" | "`svy'" == "IA42" | "`svy'" == "IA23" {
		use       `temp2', clear
		sort       Region
		local      K          = Region[_N]
		dis      "`K'"

		forvalues k = 1(1)`K' {
			dis      "`k'"
			use       `temp2', clear
			keep if    Region    == `k'
			if _N > 0 {
				local      ks         = substr("0" + string(`k'),-2,.)
				replace    fILe       = fILe + "-R" + "`ks'"
				
				drop       woman
				generate   w              = 1                                               if k == 1
				recode     w           (. = 0)
				sort       stratum cluster caseid k 
				generate   woman          = sum(w)
				drop       w
				
				keep       fILe survey country Region Urban state stratum cluster woman mother k K age sex interview B_* D_* W bidx caseid DOB Alive Born
				order      fILe survey country Region Urban state stratum cluster woman mother k K age sex interview B_* D_* W bidx caseid DOB Alive Born
				sort       fILe stratum cluster woman k
				
				export     delimited using "`pATh'/OuTPuT/`svy'-R`ks'.csv", replace		
				generate   sTArt      = year(interview)
				generate   eNd        = year(interview)
				generate   N          = 1
				collapse  (sum) N (min) sTArt (max) eNd, by(fILe survey country state)
				generate   ReTRo      = 10
				append     using `lISt'
				sort       country sTArt survey
				save      `lISt', replace
				}
			}
		}
	else {
		generate   sTArt      = year(interview)
		generate   eNd        = year(interview)	
		generate   N          = 1                                                   if mothe == 1  
		collapse  (sum) N (min) sTArt (max) eNd, by(fILe survey country)
		generate   ReTRo      = 5
		append     using `lISt'
		sort       country sTArt survey
		save      `lISt', replace	
		}
	}

merge m:1  country using "`pATh'/_codesDHS.dta", nogenerate noreport keep(master match)
sort       fILe
generate   j            = _n
order      j
export     delimited using "`pATh'/OuTPuT/lISt.csv", replace


