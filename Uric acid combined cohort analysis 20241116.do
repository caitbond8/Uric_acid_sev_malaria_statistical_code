
// Uric acid anaylsis on combined dataset 

	use "Combined_NDI_CMR01_final.dta", clear 
	
// Fig 1a
	
	bysort cohort: tab sex 
	
	tab group_cmr01 if cohort==1
	tab group_ndi  if cohort==2 
	tab smgroup cohort 
	
	tabstat age, by(cohort) stats (median)
	tab sex cohort, col
	
*------------ Only need those with uric acid data from here on ------------* 
drop if uric_se1==.

// Fig 1b
	bysort cohort: tab smgroup if uric_se1!=.
	tab smgroup if uric_ur1!=.
	
	bysort cohort: tab outcome12 smgroup
	tab follow12 smgroup 
	
	
	{ // NDI flowchart details 
	preserve 
	keep if cohort==2 
	
		tab smgroup if uric_se2!=. 
		tab smgroup if uric_se2==. 
		tab outcome if uric_se2==. 
		tab adm_death_hrs if outcome==1 & uric_se2==. 

		count if uric_se3==. & outcome==1 

		count if uric_se2==. & uric_se3!=. & smgroup==1 
		count if uric_se3==. & smgroup==1 
	
		* When children died
		count if uric_se1!=. & outcome==1 & adm_death_hrs<24  
		count if uric_se1!=. & outcome==1 & adm_death_hrs>=24 
		count if uric_se1!=. & uric_se3==. & death12m==1 & smgroup==1 

		count if uric_se1!=. & uric_se3==. & death12m!=1 & smgroup==1 & outcome!=1 

		tab follow12 if uric_se1!=. & smgroup==1 
		tab follow12 if uric_se1!=. & smgroup==0
	
	restore 
	}
	
// Fig 2a and Supplementary Table 6
	
	*Save data in excel for creating plots in Prism:
	preserve
		foreach g in 1 0 {
			foreach c in 1 2 {
				sort uric_se1 
				mkmat uric_se1 if cohort==`c' & smgroup==`g', matrix(grp`g'_`c')
			}
		}
		clear
		svmat grp1_1, name(C1_SM)
		svmat grp0_1, name(C1_CC)
		svmat grp1_2, name(C2_SM)
		svmat grp0_2, name(C2_CC)
		
		format C* %9.2f
		export excel using "prism_ready.xlsx", sheet("fig2a") ///
			firstrow(variables) sheetreplace nolabel 
	restore

	* Cross-checking plots in Stata: 
	dotplot uric_se1, over(smgroup) median bar yline(7) center, if cohort==1
	dotplot uric_se1, over(smgroup) median bar yline(7) center, if cohort==2
	
	* uric acid levels in SM vs. CC 
	bysort cohort: tabstat uric_se1, by(smgroup) stats (median p25 p75 n)
	ranksum uric_se1, by(smgroup), if cohort==1
		di `r(p)'
	ranksum uric_se1, by(smgroup), if cohort==2
		di `r(p)'
	
	* percent with elevated uric acid in SM vs CC 
	tab high_ua7 smgroup, chi exact col, if cohort==1
		di `r(p)'
	tab high_ua7 smgroup, chi exact col, if cohort==2
		di `r(p)'
	
	* SM only, over time (24 hours, 1 month)
	tabstat uric_se2 uric_se3, stats (median p25 p75 n)

	* Percent w/ elevated uric acid over time in SM
	tab high_ua7_se2 smgroup, chi exact col // % elevated UA at 24 h
	tab high_ua7_se3 smgroup, chi exact col // % elevated UA at 1 mo

	* Wilcoxon signed-rank test to compare levels over time in SM 
	signrank uric_se1 = uric_se2 // 0 to 24 h
		di `r(p)'
	signrank uric_se2 = uric_se3 // 24 h to 1 mo.
		di `r(p)'
	signrank uric_se1 = uric_se3 // 0 h to 1 mo. 
		di `r(p)'

	* McNemars test to compare proportion w/ elevated UA at 0 vs. 24 h in SM 
	mcc high_ua7 high_ua7_se2
		di `r(p)'
	mcc high_ua7_se2 high_ua7_se3
		di `r(p)'
	mcc high_ua7 high_ua7_se3
		di `r(p)'
	
// Fig 2c (Cohort 2 only) and Supplementary Table 7
	
	* Xanthine oxidase levels (serum) in SM vs CC
	tabstat xo_se1, by(smgroup) stats (median p25 p75 mean sd n)
	ranksum xo_se1, by(smgroup) 
		di `r(p)'
	  
	*Save data in excel for creating plots in Prism:
	preserve
		gen SM_xo = xo_se1_log if smgroup==1
		gen SM_uric = uric_se1_log if smgroup==1
		gen CC_xo = xo_se1_log if smgroup==0
		gen CC_uric = uric_se1_log if smgroup==0
		
		keep SM_xo SM_uric CC_xo CC_uric
		
		drop if missing(SM_xo) & missing(SM_uric) & missing(CC_xo) & missing(CC_uric)
		
		format SM* CC* %9.2f
		
		export excel using "prism_ready.xlsx", sheet("fig2c") ///
			firstrow(variables) sheetreplace nolabel
	restore
	
	* Cross-checking in Stata: 
	dotplot xo_se1_log, over(smgroup) median bar center

	*squared for quadratic regression
	foreach num of numlist 1 2 3 {
		gen uric_se`num'_sq = uric_se`num'*uric_se`num'
	}
	
	foreach num of numlist 1 2 3 {
		regress xo_se1_log  uric_se`num' uric_se`num'_sq if smgroup==1
	}

	spearman xo_se1 uric_se1 if smgroup==1
	spearman xo_se1 uric_se2 if smgroup==1
	spearman xo_se1 uric_se1 if smgroup==0
	spearman xo_se1 uric_se3 if smgroup==1
	
// Fig 2d (Cohort 2 only) and Supplementary Table 8

		
	{ // Save results table in excel 
		frame create results
		frame change results
		clear

		* Initialize empty dataset with variable names
		set obs 0
		gen str32 variable = ""
		gen double xo_rho = .
		gen double xo_pvalue = .
		gen double xo_n = .
		gen double uric_rho = .
		gen double uric_pvalue = .
		gen double uric_n = .

		frame change default

		local row = 0
		foreach var of varlist ldh_se1 hb_se1 hm_se1 cd163_pl1 hapt_pl1 hpxn_se1 hgbbase ///
			finden0 totbiomass {
			local ++row
			
			qui spearman xo_se1 `var' if smgroup==1 & cohort==2, stats (obs rho p)
			local xo_rho = r(rho)
			local xo_p = r(p)
			local xo_n = r(N)
			
			qui spearman uric_se1 `var' if smgroup==1 & cohort==2, stats (obs rho p)
			local uric_rho = r(rho)
			local uric_p = r(p)
			local uric_n = r(N)
			
			frame results {
				set obs `row'
				replace variable = "`var'" in `row'
				replace xo_rho = `xo_rho' in `row'
				replace xo_pvalue = `xo_p' in `row'
				replace xo_n = `xo_n' in `row'
				replace uric_rho = `uric_rho' in `row'
				replace uric_pvalue = `uric_p' in `row'
				replace uric_n = `uric_n' in `row'
			}
		}

		frame change results

		export excel using "Fig2_results.xlsx", sheet("2d") ///
			firstrow(variables) sheetreplace

		frame change default
		frame drop results
		}
	
// Fig 2e and Supplementary Table 9

	*Save data in excel for creating plots in Prism:
	preserve
		keep if smgroup==1 & uric_se1!=. & feua!=.
		keep uric_se1 feua
		format uric_se1 feua %9.1f
		export excel using "prism_ready.xlsx", sheet("Fig2e_1") ///
			firstrow(variables) sheetreplace nolabel
	restore

	* Stata graph cross-check	
	scatter feua uric_se1 if smgroup==1 || qfit feua uric_se1 if smgroup==1 ///
    || lowess feua uric_se1 if smgroup==1, legend(pos(6)) lcolor(red blue)
	
	*FEUA by NGAL level:
	tab ngal_ur1_150 feua_abn if smgroup==1, chi exact col

	tabstat feua, by(ngal_ur1_150) stats (median p25 p75 n), if smgroup==1
	ranksum feua, by(ngal_ur1_150), if smgroup==1
	
	*Save data in excel for creating plots in Prism:
	preserve
		sort feua
		mkmat feua if smgroup==1 & ngal_ur1_150==1, matrix(high)
		mkmat feua if smgroup==1 & ngal_ur1_150==0, matrix(low)
		clear
		svmat high, name(high_ngal)
		svmat low, name(low_ngal)		
		format high_ngal low_ngal %9.2f		
		export excel using "prism_ready.xlsx", sheet("Fig2e_2") ///
			firstrow(variables) sheetreplace nolabel
	restore

	* Cross-checking plots in Stata: 
	dotplot feua, over(ngal_ur1_150) median bar center, if smgroup==1


// Fig 2f and Supplementary Table 10

	foreach var of varlist icam_pl1 vcam_pl1 psel_pl1 esel_pl1 ang1_pl1 ang2_pl1 ///
	crp_pl1 tnfa_se1 il10_se1_pgml {
		gen `var'log=log10(`var')
	}
	
	{ // Save results table in excel 
		frame create results
		frame change results
		clear

		set obs 0
		gen str32 variable = ""
		gen cohort = .
		gen normal_mean = .
		gen normal_sd = .
		gen hyper_mean = .
		gen hyper_sd = .
		gen fold_change = .
		gen t_stat = .
		gen pvalue = .

		frame change default

		local row = 0
		foreach var of varlist icam_pl1 vcam_pl1 psel_pl1 esel_pl1 ang1_pl1 ang2_pl1 ///
			crp_pl1 tnfa_se1 il10_se1_pgml {
			foreach num of numlist 1 2 {
				qui ttest `var'log if smgroup==1 & cohort==`num', by(high_ua7)
				
				local ++row
				
				frame results {
					set obs `row'
					replace variable = "`var'" in `row'
					replace cohort = `num' in `row'
					replace normal_mean = 10^(r(mu_1)) in `row'
					replace normal_sd = 10^(r(sd_1)) in `row'
					replace hyper_mean = 10^(r(mu_2)) in `row'
					replace hyper_sd = 10^(r(sd_2)) in `row'
					replace fold_change = exp(r(mu_2)-r(mu_1)) in `row'
					replace t_stat = r(t) in `row'
					replace pvalue = r(p) in `row'
				}
			}
		}

		frame change results
		format normal_mean normal_sd hyper_mean hyper_sd fold_change %9.2f
		format t_stat %9.2f
		format pvalue %9.3f

		export excel using "Fig2_results.xlsx", sheet("2f") ///
			firstrow(variables) sheetreplace

		frame change default
		frame drop results
	}

// Supp Table 1
	foreach var of varlist bs_pos rdtresult {
		bysort cohort: tab `var' if smgroup==0 
		bysort cohort: ttest uric_se1 if smgroup==0, by(`var')
	}

// Supp Table 2


	foreach var of varlist lactate0hr ph ldh_se1 ho1_pl1 hapt_pl1 tbili_se1 ///
	hgbbase hm_se1 cd163_pl1 hpxn_se1 bun_se1 crea_se1 potass cysc_pl1 ///
	ngal_pl1 chi3l1_pl1 umod_se1 ifabp_se1 tff3_pl1 wbc plt neuabs {
		gen `var'log=log10(`var')
	}

	* Adding 1 so that values of 0 are included 
	foreach var of varlist finden0 totbiomass {
		gen `var'log = log10(`var' + 1)
	}

	{ // Save results table in excel 
		frame create results
		frame change results
		clear

		set obs 0
		gen str32 variable = ""
		gen cohort = .
		gen str32 normal_mean = ""  
		gen str32 normal_sd = ""    
		gen str32 hyper_mean = ""   
		gen str32 hyper_sd = ""     
		gen str32 fold_change = ""  
		gen str32 pvalue = ""       
		frame change default

		local row = 0
		foreach var of varlist finden0 totbiomass lactate0hr ph ldh_se1 ho1_pl1 ///
			hapt_pl1 tbili_se1 hgbbase hm_se1 cd163_pl1 hpxn_se1 bun_se1 crea_se1 potass ///
			cysc_pl1 ngal_pl1 chi3l1_pl1 umod_se1 ifabp_se1 tff3_pl1 wbc plt neuabs {
			foreach num of numlist 1 2 {
				qui count if !missing(`var'log) & smgroup==1 & cohort==`num'
				if r(N) > 0 {  // Only proceed if there are non-missing values
					qui ttest `var'log if smgroup==1 & cohort==`num', by(high_ua7)
					
					local ++row
					
					frame results {
						set obs `row'
						replace variable = "`var'" in `row'
						replace cohort = `num' in `row'
						replace normal_mean = string(10^(r(mu_1)), "%9.1f") in `row'
						replace normal_sd = string(10^(r(sd_1)), "%9.1f") in `row'
						replace hyper_mean = string(10^(r(mu_2)), "%9.1f") in `row'
						replace hyper_sd = string(10^(r(sd_2)), "%9.1f") in `row'
						replace fold_change = string(exp(r(mu_2)-r(mu_1)), "%9.2f") in `row'
						replace pvalue = string(r(p), "%9.3f") in `row'
					}
				}
				else {  // If variable is missing for this cohort
					local ++row
					
					frame results {
						set obs `row'
						replace variable = "`var'" in `row'
						replace cohort = `num' in `row'
					}
				}
			}
		}

		frame change results

		export excel using "SuppTab2.xlsx", sheet("results") ///
			firstrow(variables) sheetreplace

		frame change default
		frame drop results
		}
			
// Supp Table 3 (Cohort 2 only) 

	*Severe AKI vs. no severe AKI 
	bysort sevaki0: tabstat uric_se1 uric_se2 if smgroup==1 & uric_se2!=. , ///
		 stats (median p25 p75 n) 
	bysort sevaki0: signrank uric_se1 = uric_se2 if smgroup==1
	bysort sevaki0: tab high_ua7 if smgroup==1 & uric_se2!=. 
	bysort sevaki0: tab high_ua7_se2 if smgroup==1 

	mcc high_ua7 high_ua7_se2 if smgroup==1 & sevaki0==1
	mcc high_ua7 high_ua7_se2 if smgroup==1 & sevaki0==0
	
	gen el_ua7 = 0 if high_ua7_se2==0 & high_ua7==0
	recode el_ua7 .=1 if high_ua7_se2==0 & high_ua7==1
	recode el_ua7 .=2 if high_ua7_se2==1 & high_ua7==0
	recode el_ua7 .=3 if high_ua7_se2==1 & high_ua7==1
		label define el_ua7x 0"Not elevated at either time point" ///
		1"Elevated at admission but not 24h" ///
		2"Elevated at 24h but not at admission" ///
		3"Elevated at both time points"
		label values el_ua7 el_ua7x

	tab el_ua7 sevaki0, chi exact col

	*AKI vs no AKI
	bysort aki0: tabstat uric_se1 uric_se2 if smgroup==1 &  uric_se2!=. , ///
		stats (median p25 p75 n) 
	bysort aki0: signrank uric_se1 = uric_se2 if smgroup==1
	bysort aki0: tab high_ua7 if smgroup==1 & uric_se2!=. 
	bysort aki0: tab high_ua7_se2 if smgroup==1 & uric_se1!=. 

	mcc high_ua7 high_ua7_se2 if smgroup==1 & aki0==1
	mcc high_ua7 high_ua7_se2 if smgroup==1 & aki0==0
	
	tab el_ua7 aki0, chi exact col

	
// Extended Data Fig 1d

	* compare uric acid levels by genotype stratified by SM group 
	graph box uric_se1, over(chr49942428AG_A) by(smgroup) yline(7) legend(col(1)) 
	graph box uric_se3, over(chr49942428AG_A)yline(7) legend(col(1)) 

	ttest uric_se1, by(chr49942428_G), if smgroup==1
		
	tabstat uric_se1, by(chr49942428AG_A) stats (median p25 p75 n), if smgroup==1
	kwallis uric_se1, by(chr49942428AG_A), if smgroup==1
	dunntest uric_se1, by(chr49942428AG_A), if smgroup==1
	nptrend uric_se1, group(chr49942428AG_A) cuzick , if smgroup==1
	
	tabstat uric_se3, by(chr49942428AG_A) stats (median p25 p75 n), if smgroup==1
	kwallis uric_se3, by(chr49942428AG_A), if smgroup==1
	nptrend uric_se3, group(chr49942428AG_A) cuzick , if smgroup==1
		
// Supp Table 4a
	
	tab chr49942428_G outcome, chi exact col 
	tab chr49942428AG_A outcome, chi exact col 

// Supp Table 4b

	logistic outcome uric_se1 if chr49942428AG_A==0
	logistic outcome uric_se1 if chr49942428AG_A==1
	logistic outcome uric_se1 if chr49942428AG_A==2

// Extended Data Fig 3b 
	
	tab int_inj high_ua7 , chi exact col, if smgroup==1
	
	* See R file "Uric_Microbiome" for Fig 3c and 3d analysis

// Fig 3e and Supp Table 19	

	* Association between bacteria (Relative abundance) and intestinal injury markers
	
	{ // Save results table in excel 
	frame create results
	frame change results
	clear

	* Initialize empty dataset with variable names
	set obs 0
	gen str32 variable = ""
	gen str32 no_injury = ""  
	gen str32 injury = ""      
	gen z_stat = .
	gen pvalue = .

	frame change default

	* Initialize row counter
	local row = 0

	* Loop through variables in the specified order
	foreach var of varlist escherichia shigella enterobacter klebsiella enterococcus {
		local ++row
		
		* Get statistics for no injury (int_inj=0)
		qui sum `var' if int_inj==0 & smgroup==1, detail
		local med0 = string(r(p50), "%3.1f")
		local p25_0 = string(r(p25), "%3.1f")
		local p75_0 = string(r(p75), "%3.1f")
		
		* Get statistics for injury (int_inj=1)
		qui sum `var' if int_inj==1 & smgroup==1, detail
		local med1 = string(r(p50), "%3.1f")
		local p25_1 = string(r(p25), "%3.1f")
		local p75_1 = string(r(p75), "%3.1f")
		
		* Get ranksum statistics
		qui ranksum `var' if smgroup==1, by(int_inj)
		
		frame results {
			set obs `row'
			replace variable = "`var'" in `row'
			replace no_injury = "`med0' (`p25_0', `p75_0')" in `row'
			replace injury = "`med1' (`p25_1', `p75_1')" in `row'
			replace z_stat = r(z) in `row'
			replace pvalue = r(p) in `row'
		}
	}

	frame change results

	format z_stat %9.3f
	format pvalue %9.3f

	export excel using "ED_Fig2e_results.xlsx", sheet("results") ///
		firstrow(variables) sheetreplace

	frame change default
	frame drop results

	}
		
	* Save data to excel for Prism figures
	preserve
		sort escherichia
		mkmat escherichia if int_inj==1 & smgroup==1, matrix(esc1)
		mkmat escherichia if int_inj==0 & smgroup==1, matrix(esc0)
		
		sort shigella
		mkmat shigella if int_inj==1 & smgroup==1, matrix(shi1)
		mkmat shigella if int_inj==0 & smgroup==1, matrix(shi0)
		
		sort klebsiella
		mkmat klebsiella if int_inj==1 & smgroup==1, matrix(kleb1)
		mkmat klebsiella if int_inj==0 & smgroup==1, matrix(kleb0)
		
		sort enterobacter
		mkmat enterobacter if int_inj==1 & smgroup==1, matrix(ent1)
		mkmat enterobacter if int_inj==0 & smgroup==1, matrix(ent0)
		
		sort bacteroides
		mkmat bacteroides if int_inj==1 & smgroup==1, matrix(bac1)
		mkmat bacteroides if int_inj==0 & smgroup==1, matrix(bac0)
		
		sort enterococcus
		mkmat enterococcus if int_inj==1 & smgroup==1, matrix(enc1)
		mkmat enterococcus if int_inj==0 & smgroup==1, matrix(enc0)
		
		clear
		svmat esc1, name(Escherichia_IntInj)
		svmat esc0, name(Escherichia_NoIntInj)
		svmat shi1, name(Shigella_IntInj)
		svmat shi0, name(Shigella_NoIntInj)
		svmat kleb1, name(Klebsiella_IntInj)
		svmat kleb0, name(Klebsiella_NoIntInj)
		svmat ent1, name(Enterobacter_IntInj)
		svmat ent0, name(Enterobacter_NoIntInj)
		svmat bac1, name(Bacteroides_IntInj)
		svmat bac0, name(Bacteroides_NoIntInj)
		svmat enc1, name(Enterococcus_IntInj)
		svmat enc0, name(Enterococcus_NoIntInj)
		
		export excel using "prism_ready.xlsx", sheet("ED_Fig2e") ///
			firstrow(variables) sheetreplace

	restore

	
// Fig 3a and Supp Table 11

	* For Prism: 
	preserve
		* Create matrices for each group
		sort uric_se1
		mkmat uric_se1 if smgroup==1 & outcome==1 & cohort==1, matrix(death_c1)
		mkmat uric_se1 if smgroup==1 & outcome==0 & cohort==1, matrix(alive_c1)
		mkmat uric_se1 if smgroup==1 & outcome==1 & cohort==2, matrix(death_c2)
		mkmat uric_se1 if smgroup==1 & outcome==0 & cohort==2, matrix(alive_c2)
		
		* Create new dataset with these columns
		clear
		svmat death_c1, name(Cohort1_Death)
		svmat alive_c1, name(Cohort1_Alive)
		svmat death_c2, name(Cohort2_Death)
		svmat alive_c2, name(Cohort2_Alive)
		
		format Cohort* %9.1f
		
		* Export to Excel
		export excel using "prism_ready.xlsx", sheet("Fig3a") ///
			firstrow(variables) sheetreplace
	restore

	* Stata cross-check 
	dotplot uric_se1, over(outcome) median bar yline(7) center, if smgroup==1 & cohort==1
	dotplot uric_se1, over(outcome) median bar yline(7) center, if smgroup==1 & cohort==2 

	* Uric acid elvels in those who died vs. survived
	bysort cohort: tabstat uric_se1, by(outcome) stats (median p25 p75 n), if smgroup==1
	ranksum uric_se1, by(outcome),  if smgroup==1 & cohort==1
		di `r(p)'
	ranksum uric_se1, by(outcome),  if smgroup==1 & cohort==2
		di `r(p)'

	* Percent w/ elevated UA in those who died vs. survived 
	tab high_ua7 outcome, chi exact col, if smgroup==1 & cohort==1
		di `r(p)'
	tab high_ua7 outcome, chi exact col, if smgroup==1 & cohort==2
		di `r(p)'

	* Uric acid levels at 24h in those who died vs. survived
	tabstat uric_se2, by(outcome) stats (median p25 p75 n), if smgroup==1
	ranksum uric_se2, by(outcome),  if smgroup==1
		di `r(p)'

	* Percent w/ elevated UA in those who died vs. survived 
	tab high_ua7_se2 outcome, chi exact col, if smgroup==1 
		di `r(p)'

// Fig 3b and Supp Table 12

	gen uric_se1_ln = log(uric_se1)	
	gen uric_se2_ln = log(uric_se2)	

cls
	foreach var of varlist high_ua7 uric_se1_ln {
		
		** Crude
		* Cohort 1
		logistic outcome `var' if smgroup==1 & cohort==1
		test `var'
		display "Cohort 1 - `var' - `r(p)'"
		* Cohort 2
		logistic outcome `var' if smgroup==1 & cohort==2
		test `var'
		display "Cohort 2 - `var' - `r(p)'"
		
		** Adjusted
		*Cohort 1
		logistic outcome `var' age sex if smgroup==1 & cohort==1 
		test `var'
		display "Cohort 2 - `var' - `r(p)'"
		*Cohort 2 (also adjusting for site)
		logistic outcome `var' age sex site if smgroup==1 & cohort==2
		test `var'
		display "Cohort 2 - `var' - `r(p)'"
	}
	
	
cls
	foreach var of varlist high_ua7_se2 uric_se2_ln {
		** Crude
		logistic outcome `var' if smgroup==1 & cohort==2
		test `var'
		display "Cohort 2 - `var' - `r(p)'"
		
		** Adjusted
		logistic outcome `var' age sex site if smgroup==1 & cohort==2
		test `var'
		display "Cohort 2 - `var' - `r(p)'"
		}
		
// Fig 3c and Supp Table 13

	cls
	foreach var of varlist lacidosis aki0 coma  {
		logistic high_ua7 `var' if smgroup==1 & cohort==1
		test `var'
		display "Cohort 1 - `var' - `r(p)'"
	}		

	foreach var of varlist lacidosis aki0 coma acidosis_ndi coldperi int_inj {

		logistic high_ua7 `var' if smgroup==1 & cohort==2
		test `var'
		display "Cohort 2 - `var' - `r(p)'"
	}			
	
	
	
	* Cohort 1 - Significant from univariable model: lacidosis and aki
	logistic high_ua7 lacidosis aki0 age sex if smgroup==1 & cohort==1 
	test aki0 
	di `r(p)' // aki

	* Cohort 2 - Significant from univariable model: aki, coma, acidosis, coldperi, int_inj
	logistic high_ua7 aki0 coma acidosis_ndi coldperi int_inj age sex site if smgroup==1 & cohort==2
	test aki0 
	di `r(p)' // aki
	test acidosis_ndi 
	di `r(p)' // acidosis 
	test int_inj 
	di `r(p)' // intestinal injury


// Fig 3d (Cohort 2 only) 

	foreach var of varlist lacidosis aki0 coma acidosis coldperi int_inj {
		tab high_ua7 `var', chi exact col, if smgroup==1 & cohort==2 
	}	
		
// Fig 3e and Supp Table 14 (Cohort 2 only) 

	{ // Excel output for proportions and 95% CIs: 
	frame create results
	frame change results
	clear
	set obs 0
	gen str32 condition = ""
	gen str32 category = ""
	gen str32 probability = ""
	frame change default
	local row = 0
	foreach var of varlist lacidosis aki0 coma acidosis coldperi int_inj {
		qui logistic outcome `var'##high_ua7 age sex site if cohort==2 
		qui margins `var'#high_ua7
		matrix m = r(b)
		matrix v = r(V)
		local g1 "No Normal uric acid"
		local g2 "No Hyperuricemia"
		local g3 "Yes Normal uric acid"
		local g4 "Yes Hyperuricemia"	
		forvalues i = 1/4 {
			local ++row
			local prob = 100*m[1,`i']
			local se = sqrt(v[`i',`i'])
			local lb = 100*(m[1,`i'] - 1.96*`se')
			local ub = 100*(m[1,`i'] + 1.96*`se')
			frame results {
				set obs `row'
				replace condition = "`var'" in `row'
				replace category = "`g`i''" in `row'
				replace probability = string(`prob', "%3.1f") + "% (" + string(`lb', "%3.1f") + ", " + string(`ub', "%3.1f") + ")" in `row'			
			}
		}	
		local ++row
		frame results {
			set obs `row'
		}
	}
	frame change results
	export excel using "Fig3e_results.xlsx", sheet("prob") ///
		firstrow(variables) sheetreplace
	}
			
			
	{ // Excel output for Sidak pairwise p-values: 
	frame change default
	frame drop results
	frame create results
	frame change results
	clear
	set obs 0
	gen str32 condition = ""
	gen double sidak1 = .
	gen double sidak2 = .
	gen double sidak3 = .
	gen double sidak4 = .
	gen double sidak5 = .
	gen double sidak6 = .
	frame change default
	local row = 0
	foreach var of varlist lacidosis aki0 coma acidosis_ndi coldperi int_inj {
		local ++row
		qui logistic outcome `var'##high_ua7 age sex site if cohort==2 
		qui pwcompare `var'#high_ua7, effects sid
		matrix p = r(table_vs)
		frame post results ("`var'") ///
			(p[4,1]) /// No#Elevated UA vs No#Normal UA
			(p[4,2]) /// Yes#Normal UA vs No#Normal UA
			(p[4,4]) /// Yes#Elevated UA vs No#Normal UA
			(p[4,5]) /// Yes#Normal UA vs No#Elevated UA
			(p[4,3]) /// Yes#Elevated UA vs No#Elevated UA
			(p[4,6])  // Yes#Elevated UA vs Yes#Normal UA
	}
	frame change results
	format sidak* %9.3f
	export excel using "Fig3e_results.xlsx", sheet("sidak_pvalues") ///
		firstrow(variables) sheetreplace
	frame change default
	frame drop results
	
	* Add a note about sidak p values: 
	putexcel set "Fig3e_results.xlsx", modify sheet("sidak_pvalues")
		putexcel A10 = "Comparison numbers corresponding with sidak p-values"
		putexcel A11 = "1: No/HU vs. No/normal"
		putexcel A12 = "2: Yes/normal vs. No/normal"
		putexcel A13 = "3: Yes/HU vs. No/normal" 
		putexcel A14 = "4: Yes/normal vs. No/HU"
		putexcel A15 = "5: Yes/HU vs. No/HU"
		putexcel A16 = "6: Yes/HU vs. Yes/normal"
		putexcel A18 = "HU = hyperuricemia"
	}
	
	cls
	* Numbers died / total N per group (for figure Fig3e)
	foreach var of varlist lacidosis aki0 coma acidosis coldperi int_inj {
		di ""
		di "Neither `var' or hyperuricemia"
		tab outcome if `var'==0 & high_ua7==0 & cohort==2 
		di ""
		di "Hyperuricemia but no `var'"
		tab outcome if `var'==0 & high_ua7==1 & cohort==2   
		di ""
		di "`var' but no hyperuricemia"
		tab outcome if `var'==1 & high_ua7==0 & cohort==2 
		di ""
		di "Both `var' and hyperuricemia"
		tab outcome if `var'==1 & high_ua7==1 & cohort==2  
	}
	

// Extended data Fig 3 (Cohort 2 only) 

	lowess outcome uric_se1 if cohort==2 
	
	*Make splines based on change in direction (lowess) and hyperuricemia cutoff 7	
	mkspline ua1 3.8 ua2 7 ua3 = uric_se1 if cohort==2
	logistic outcome ua1 ua2 ua3 age sex site if cohort==2
	
// Extended data Fig 4 (Cohort 2 only) 

cls
preserve 
keep if cohort==2 

	foreach m of varlist aki0 coma acidosis coldperi int_inj {
		di""
		di"**** `m'****" 
		mediate (outcome age sex site, logit) (`m' age sex site, logit) (high_ua7), nointeraction
		di ""
		di "****************************"
		di "Percent mediated = " = (r(table)[1,1]/r(table)[1,3])*100
		di "****************************"		
		logistic outcome `m' high_ua7 age sex site
		logistic `m' high_ua7 age sex site if smgroup==1
	}
		
	foreach x of varlist aki0 { 
		di""
		di"**** `x'****" 
		mediate (outcome age sex site, logit) (high_ua7 age sex site, logit) (`x'), nointeraction
		di ""
		di "****************************"
		di "Percent mediated = " = (r(table)[1,1]/r(table)[1,3])*100
		di "****************************"		
		logistic outcome `x' high_ua7 age sex site
		logistic high_ua7 `x' age sex site if smgroup==1
	}	
	
restore 
	
	
// Fig 4a and Supp Table 15 

	*Save data in excel for creating plots in Prism:
	preserve
	keep if smgroup==1 & outcome!=1
	
		foreach g in 1 0 {
			foreach c in 1 2 {
				sort uric_se1 
				mkmat uric_se1 if cohort==`c' & death12m==`g', matrix(grp`g'_`c')
			}
		}
		clear
		svmat grp1_1, name(C1_Died)
		svmat grp0_1, name(C1_Survived)
		svmat grp1_2, name(C2_Died)
		svmat grp0_2, name(C2_Survived)
		
		format C* %9.2f
		export excel using "prism_ready.xlsx", sheet("fig4a_0h") ///
			firstrow(variables) sheetreplace nolabel 
	restore
	
	*24 hours uric acid (only in Cohort 2)
	preserve
	keep if smgroup==1 & outcome!=1
	
		foreach g in 1 0 {
			foreach c in 1 2 {
				sort uric_se2 
				mkmat uric_se2 if cohort==`c' & death12m==`g', matrix(grp`g'_`c')
			}
		}
		clear
		svmat grp1_1, name(C1_Died)
		svmat grp0_1, name(C1_Survived)
		svmat grp1_2, name(C2_Died)
		svmat grp0_2, name(C2_Survived)
		
		format C* %9.2f
		export excel using "prism_ready.xlsx", sheet("fig4a_24h") ///
			firstrow(variables) sheetreplace nolabel 
	restore
	
	
	*gen uric_se1_ln=log(uric_se1) 

	cls
	preserve 
	keep if smgroup==1
		*0h
		bysort cohort: tabstat uric_se1, by(death12m) stats (median p25 p75 n)
		ranksum uric_se1, by(death12m) , if cohort==1
			di `r(p)'
		ranksum uric_se1, by(death12m) , if cohort==2
			di `r(p)'

		tab high_ua7 death12m, chi exact col, if cohort==1
			di `r(p)'
		tab high_ua7 death12m, chi exact col, if cohort==2
			di `r(p)'
		
		*24h
		tabstat uric_se2, by(death12m) stats (median p25 p75 n)
		ranksum uric_se2, by(death12m)

		tab high_ua7_se2 death12m, chi exact col
		
	restore
		
// Fig 4b and Supp Table 16

	preserve 
	keep if smgroup==1

	* Cohort 1 (use Firth logistic regression since sex perfectly predicts outcome)
	foreach var of varlist uric_se1_ln  high_ua7  {
		firthlogit death12m `var' if cohort==1, or
		firthlogit death12m `var' age sex  if cohort==1, or		
	}	
	* Cohort 2
	foreach var of varlist uric_se1_ln uric_se2_ln high_ua7 high_ua7_se2 {
		
		logistic death12m `var' if cohort==2
		logistic death12m `var' age sex site if cohort==2
	}
		
	restore 


// Fig 4c
	
	stset deathdaysto12m_all if smgroup==1 & outcome!=1, failure(death12m)
	bysort cohort: stcox high_ua7 if smgroup==1

	*Cohort 1
	sts graph, by(high_ua7) risktable  tmax(354) xlabel(0 (73) 365), ///
	if smgroup==1 & outcome!=1 & cohort==1 	
	sts, by(high_ua7)   yla(.75 .8 .85 .9 .95 1, angle(0)) tmax(372) ///
	xlabel(0 (73) 372) , if smgroup==1 & outcome!=1 & cohort==1 
	sts test high_ua7 if cohort==1 // log-rank test 

	*Cohort 2
	sts graph, by(high_ua7) risktable  tmax(354) xlabel(0 (73) 365), ///
	if smgroup==1 & outcome!=1 & cohort==2
	sts, by(high_ua7)   yla(.75 .8 .85 .9 .95 1, angle(0)) tmax(372) ///
	xlabel(0 (73) 372) , if smgroup==1 & outcome!=1 & cohort==2 
	sts test high_ua7 if cohort==2 // log-rank test 

	*For Prism - Cohort 1
	preserve
		keep if smgroup==1 & outcome==0 & cohort==1
		keep deathdaysto12m_all death12m high_ua7  
		gen death_normal = death12m if high_ua7==0
		gen death_hyper = death12m if high_ua7==1
		keep deathdaysto12m_all death_normal death_hyper
		export excel using "prism_ready.xlsx", sheet("Fig4c_c1") ///
			firstrow(variables) sheetreplace
	restore

	*For Prism - Cohort 2
	preserve
		keep if smgroup==1 & outcome!=1 & cohort==2
		keep deathdaysto12m_all death12m high_ua7
		gen death_normal = death12m if high_ua7==0
		gen death_hyper = death12m if high_ua7==1
		keep deathdaysto12m_all death_normal death_hyper
		export excel using "prism_ready.xlsx", sheet("Fig4c_c2") ///
			firstrow(variables) sheetreplace
	restore

// Fig 4d

	* Venn diagram 
	bysort cohort: count if sev_anemia==1 & high_ua7==1 & death12m!=.
	bysort cohort: count if sev_anemia==1 & high_ua7==0 & death12m!=.
	bysort cohort: count if sev_anemia==0 & high_ua7==1 & death12m!=.
	bysort cohort: count if sev_anemia==0 & high_ua7==0 & death12m!=.

	bysort cohort: tab death12 if sev_anemia==0 & high_ua7==0 & death12m!=.
	bysort cohort: tab death12 if sev_anemia==1 & high_ua7==0 & death12m!=.
	bysort cohort: tab death12 if sev_anemia==0 & high_ua7==1 & death12m!=.
	bysort cohort: tab death12 if sev_anemia==1 & high_ua7==1 & death12m!=.
	

	*Cohort 1
	logistic death12m i.high_ua7##i.sev_anemia if cohort==1
	margins high_ua7#sev_anemia 
	pwcompare high_ua7#sev_anemia, effects sid

	*Cohort 2
	logistic death12m i.high_ua7##i.sev_anemia if cohort==2
	margins high_ua7#sev_anemia 
	pwcompare high_ua7#sev_anemia, effects sid

// Extended Data Fig 5  -- Flowchart for cognitive outcomes 

cls


preserve 
keep if cohort==1  // Cohort 1 

	tab outcome if smgroup==1, mi 
 
	tab outcome if smgroup==1 & cogdat0_u5==1 
	tab outcome if smgroup==1 & cogdat0_a5==1 
	
	tab outcome6mo if smgroup==1 &  outcome==0
 
	tab outcome6mo if smgroup==1 &  outcome==0 & cogdat6_u5==1 
	tab outcome6mo if smgroup==1 &  outcome==0 & cogdat6_a5==1

	tab outcome12mo if smgroup==1 &  outcome==0 & outcome6mo==0
 
	tab outcome12mo if smgroup==1 &  outcome==0 & outcome6mo==0 & cogdat12_u5==1 
	tab outcome12mo if smgroup==1 &  outcome==0 & outcome6mo==0 & cogdat12_a5==1

	tab outcome24mo if smgroup==1 &  outcome==0 & outcome6mo==0 ///
	& outcome12mo==0
 
	tab outcome24mo if smgroup==1 &  outcome==0 & outcome6mo==0 ///
	& outcome12mo==0 & cogdat24_u5==1
	tab outcome24mo if smgroup==1 &  outcome==0 & outcome6mo==0 ///
	& outcome12mo==0 & cogdat24_a5==1 
	
restore 
	
	
	// Cohort 2
	
	tab smgroup if cohort==2
	count if smgroup==1 & cogdat12_u5==1 & cohort==2
	tab follow12 if smgroup==1  & cogdat12_u5!=1 & cohort==2
		

// Fig 5a - Cohort 1 (see below, running last due to data restructure)

// Fig 5a - Cohort 2 
	
	putexcel set Fig5, modify sheet("5a_C2")
	putexcel C1=("N") D1=("Beta, (95% CI)") E1=("P")

	local row = 2
	preserve 
	keep if cohort==2 

	foreach contvar1 of varlist mulout12_zscr ecvt12_zscr coatsc12_zscr {
		// Unadjusted
		putexcel A`row'=("`contvar1'")
		putexcel B`row'=("Unadjusted")
		
		regress `contvar1' uric_se1_log if smgroup==1
		
		mat r = r(table)
		if r[4,1] <= 0.05 putexcel E`row', fpattern(solid, lightyellow)
		
		putexcel C`row'=("`e(N)'") ///
			D`row'=("`:display %4.2f r[1,1]' (`:display %4.2f r[5,1]', `:display %4.2f r[6,1]')") ///
			E`row'=("`:display %4.3f r[4,1]'"), hcenter
		
		// Adjusted
		local ++row
		putexcel B`row'=("Adjusted")
		
		regress `contvar1' uric_se1_log age sex site cm waz0 haz0 sqrt_malariafu if smgroup==1
		
		mat r = r(table)
		if r[4,1] <= 0.05 putexcel E`row', fpattern(solid, lightyellow)
		
		putexcel C`row'=("`e(N)'") ///
			D`row'=("`:display %4.2f r[1,1]' (`:display %4.2f r[5,1]', `:display %4.2f r[6,1]')") ///
			E`row'=("`:display %4.3f r[4,1]'"), hcenter
		
		local row = `row' + 2
	}
	restore


// Fig 5b - Cohort 1 (continued in long format below)

	putexcel set Fig5, modify sheet("5b")

	// Set up headers
	putexcel A1=("Biomarker") B1=("Normal uric acid") C1=("Hyperuricemia") ///
		D1=("z") E1=("P-value")

	local row = 2
	foreach var of varlist mda_cs1 sodact_cs1 sodconc_cs1 adma_cs1 tau_cs1 kynu_cs1 kyna_cs1 {
		// Get medians and IQRs for each group
		quietly tabstat `var', by(high_ua7) stats(median p25 p75) save, if smgroup==1
		matrix stats = r(Stat1), r(Stat2)
		
		// Get ranksum test results
		quietly ranksum `var', by(high_ua7), if smgroup==1
		local pvalue = r(p)
		local z = r(z)
		
		// Write to Excel
		putexcel A`row'=("`var'") ///
			B`row'=("`:display %4.2f stats[1,1]' (`:display %4.2f stats[2,1]', `:display %4.2f stats[3,1]')") ///
			C`row'=("`:display %4.2f stats[1,2]' (`:display %4.2f stats[2,2]', `:display %4.2f stats[3,2]')") ///
			D`row'=("`:display %4.3f `z''") ///
			E`row'=("`:display %4.3g `pvalue''"), hcenter
			
		if `pvalue' <= 0.05 putexcel E`row', fpattern(solid, lightyellow)
		
		local ++row
	}

	*graph cross check
	qui {
		foreach var of varlist mda_cs1 sodact_cs1 sodconc_cs1 adma_cs1 ///
		tau_cs1	kynu_cs1 kyna_cs1   {
			graph box `var', over(high_ua7) name(g_`var', replace) title(`var') 
		}
	}
		// Combined 
		graph combine g_mda_cs1 g_sodact_cs1 g_sodconc_cs1 g_adma_cs1 ///
		g_tau_cs1 g_kynu_cs1 g_kyna_cs1 


	*for write up
	tab alb_modsev high_ua7, chi exact col

	* For prism to make box plots: 
	preserve
		sort mda_cs1
		mkmat mda_cs1 if high_ua7==1, matrix(mda1)
		mkmat mda_cs1 if high_ua7==0, matrix(mda0)
		sort sodact_cs1
		mkmat sodact_cs1 if high_ua7==1, matrix(sodact1)
		mkmat sodact_cs1 if high_ua7==0, matrix(sodact0)
		sort sodconc_cs1
		mkmat sodconc_cs1 if high_ua7==1, matrix(sodconc1)
		mkmat sodconc_cs1 if high_ua7==0, matrix(sodconc0)
		sort tau_cs1
		mkmat tau_cs1 if high_ua7==1, matrix(tau1)
		mkmat tau_cs1 if high_ua7==0, matrix(tau0)	
		sort kynu_cs1
		mkmat kynu_cs1 if high_ua7==1, matrix(kynu1)
		mkmat kynu_cs1 if high_ua7==0, matrix(kynu0)
		sort kyna_cs1
		mkmat kyna_cs1 if high_ua7==1, matrix(kyna1)
		mkmat kyna_cs1 if high_ua7==0, matrix(kyna0)	
		sort adma_cs1
		mkmat adma_cs1 if high_ua7==1, matrix(adma1)
		mkmat adma_cs1 if high_ua7==0, matrix(adma0)
		
		clear
		svmat mda1, name(MDA_Hyperuricemia)
		svmat mda0, name(MDA_Normal)
		svmat sodact1, name(SODact_Hyperuricemia)
		svmat sodact0, name(SODact_Normal)
		svmat sodconc1, name(SODconc_Hyperuricemia)
		svmat sodconc0, name(SODconc_Normal)
		svmat tau1, name(Tau_Hyperuricemia)
		svmat tau0, name(Tau_Normal)
		svmat kynu1, name(Kynu_Hyperuricemia)
		svmat kynu0, name(Kynu_Normal)
		svmat kyna1, name(Kyna_Hyperuricemia)
		svmat kyna0, name(Kyna_Normal)
		svmat adma1, name(ADMA_Hyperuricemia)
		svmat adma0, name(ADMA_Normal)
		
		export excel using "prism_ready.xlsx", sheet("Fig5b") ///
			firstrow(variables) sheetreplace
	restore
	

// Supp Table 5e

preserve 
keep if cohort==2 

	* Watson's reserach criteria only: 
	foreach contvar1 of varlist mulout12_zscr ecvt12_zscr coatsc12_zscr {
		regress `contvar1' uric_se1_log if smgroup==1  & watson_sm==1
		regress `contvar1' uric_se1_log age sex site cm waz0 haz0 sqrt_malariafu if smgroup==1  & watson_sm==1
	}
	
restore 	
	

	
	
	
********************************************************************************

// Fig 5a - Cohort 1 - run this last 
	
	// First run new z-score calculations (use full cohort for z-score calcs)
	use "Combined_NDI_CMR01_final.dta", clear 
	
	drop if cohort==2
		
		keep studyid smgroup group_cmr01 sex agefin* waz* haz* sqrt_malariafu ///
		mulout* ecvt* coatsc* mpi* vdpt* seqpro* agefin* treatgrp uric* agecogt ///
		mda_cs1  sodact_cs1 sodconc_cs1 adma_cs1 tau_cs1 kynu_cs1 kyna_cs1  ///
		watson_sm 

		reshape long agefin mulout ecvt coatsc mpi vdpt seqpro, i(studyid) j(visit)
		
		*excluding scores that do not meet criteria 
		replace mulout = . if visit==6 & (studyid==47 | studyid==24 | studyid==9 | studyid==2)
		replace mulout = . if visit==12 & (studyid==100 | studyid==2 | studyid==86 | studyid==79)
		replace ecvt = . if visit==0 & studyid==82
		replace ecvt = . if visit==6 & (studyid==47 | studyid==24 | ///
		studyid==9 | studyid==2 | studyid==319 | studyid==359)
		replace ecvt = . if visit==12 & (studyid==8 | studyid==100 | studyid==1 ///
		| studyid==2 | studyid==86 | studyid==79 | studyid==110 | studyid==47 | ///
		studyid==103 | studyid==62 | studyid==52 | studyid==9 | studyid==50 | studyid==232)
		replace coatsc = . if visit==6 & (studyid==47 | studyid==24 | studyid==9 | studyid==2)
		replace coatsc = . if visit==12 & (studyid==100 | studyid==2 | studyid==86 | studyid==79)

		* need to calculate mean age for CC (not U01) for each test 
		foreach var of varlist mulout ecvt coatsc mpi vdpt seqpro {
			sum agefin if `var'!=. & smgroup==0 & treatgrp==.
			gen age_mean_`var' = r(mean)
			gen age_cent_`var' = agefin - age_mean_`var' 
			gen age_cent2_`var' = age_cent_`var'^2
		}
	
		* For U01 calcualtions, use CC in U01 who did not receive iron (treatgrp==3)
		foreach var of varlist mulout ecvt coatsc mpi vdpt seqpro {
			sum agefin if `var'!=. & smgroup==0 & treatgrp==3
			gen age_mean_`var'_u01 = r(mean)
			gen age_cent_`var'_u01 = agefin - age_mean_`var' 
			gen age_cent2_`var'_u01 = age_cent_`var'^2
		}	
	
	// Run model to calculate z-scores 

		cls
		foreach var of varlist mulout ecvt coatsc mpi vdpt seqpro {
			mixed `var' age_cent_`var' age_cent2_`var' if smgroup==0 & `var'!=. & treatgrp==. ///
			|| studyid: , reml covariance(unstructured) dfmethod(kroger) iterate(50) 

				predict `var'_mean, xb // average score 	
				
				matrix b = e(b)
				local id_var_col = colnumb(b, "lns1_1_1:_cons")
				local id_variance = exp(b[1, `id_var_col'])^2 // random intercept var
				local res_var_col = colnumb(b, "lnsig_e:_cons")
				local res_variance = exp(b[1, `res_var_col'])^2 // residual var
				local `var'_yvar = sqrt(`id_variance' + `res_variance') // Y SD
			
				gen `var'_zscr = (`var' - `var'_mean)/ ``var'_yvar' if `var'!=. & treatgrp==.
				
				di""
				di "Variance of Y =" ``var'_yvar'
		}
		
		
		cls
		foreach var of varlist mulout ecvt coatsc mpi vdpt seqpro {
			mixed `var' age_cent_`var'_u01 age_cent2_`var'_u01 if smgroup==0 & `var'!=. & treatgrp==3 ///
			|| studyid: , reml covariance(unstructured) dfmethod(kroger) iterate(50) 
				predict `var'_mean_u01, xb // average score 	
				
				matrix b = e(b)
				local id_var_col = colnumb(b, "lns1_1_1:_cons")
				local id_variance = exp(b[1, `id_var_col'])^2 // random intercept var
				local res_var_col = colnumb(b, "lnsig_e:_cons")
				local res_variance = exp(b[1, `res_var_col'])^2 // residual var
				local `var'_yvar = sqrt(`id_variance' + `res_variance') // Y SD
			
				replace `var'_zscr = (`var' - `var'_mean_u01)/ ``var'_yvar' if `var'!=. & treatgrp!=.
				
				di""
				di "Variance of Y =" ``var'_yvar'
		}
			
		foreach var of varlist mulout ecvt coatsc mpi vdpt seqpro {
			replace `var'_zscr = -5 if `var'_zscr<-5 & `var'_zscr!=.
			replace `var'_zscr = 5 if `var'_zscr>5 & `var'_zscr!=.
		}	
		
	generate sgroup = .
	replace sgroup = 0 if group_cmr01 == 2
	replace sgroup = 1 if group_cmr01 == 1
		
	
	// Adjusted and adjusted LME
	keep if smgroup==1

	* write a program for the LME models
	capture program drop run_mixed_analysis
	capture program drop write_results
	matrix drop _all

	program define run_mixed_analysis, rclass
		args var condition adjust
		
		// Build the full model command
		local model "`var'_zscr uric_se1_log"
		if "`adjust'" == "adjusted" {
			local model "`model' agefin00 sex waz0 haz0 sqrt_malariafu sgroup"
		}
		local model "`model' i.visit"
		if "`condition'" != "" {
			local model "`model' if `condition'"
		}
		
		mixed `model' || studyid:, reml covariance(unstructured) ///
			residuals(banded 0, t(visit)) dfmethod(kroger) iterate(50)
		
		matrix values = r(table)
		local orb = values[1,1]
		local cilow = values[5,1]
		local cihigh = values[6,1]
		local adpvalue = values[4,1]
		matrix ns = e(N_g)    
		local adn = ns[1,1]
		local adnobs = e(N)
		estat ic
		matrix table = r(S)
		local aic = table[1,5]
		
		// Store results in global matrix
		matrix results = (`adnobs', `adn', `orb', `cilow', `cihigh', `adpvalue', `aic')
		matrix colnames results = n_obs n beta ci_low ci_high p_value aic
		
		display "Beta: `orb', P-value: `adpvalue'"
	end

	program define write_results
		args row var
		
		// Get values from the global matrix 'results'
		local adnobs = results[1,1]
		local adn = results[1,2]
		local orb = results[1,3]
		local cilow = results[1,4]
		local cihigh = results[1,5]
		local adpvalue = results[1,6]
		local aic = results[1,7]
		
		if `adpvalue' <= 0.05 {
			putexcel E`row', fpattern(solid, lightyellow)
		}
		
		putexcel B`row'=("`var'") ///
			C`row'=("`adnobs', `adn'") ///
			D`row'=("`:display %4.2f `orb'' (`:display %4.2f `cilow'', `:display %4.2f `cihigh'')") ///
			E`row'=("`:display %4.3f `adpvalue''") ///
			G`row'=("`:display %4.2f `orb''") ///
			H`row'=("`:display %4.2f `cilow''") ///
			I`row'=("`:display %4.2f `cihigh''") ///
			K`row'=("`:display %4.2f `aic''"), hcenter
	end

	// adjusted 
	putexcel set Fig5, modify sheet("5a_C1_Adjusted")

	local row = 1
	putexcel C`row'=("N(obs), n") D`row'=("Beta, (95% CI)") E`row'=("P")
	putexcel F`row'=("Adjusted") G`row'=("Beta") H`row'=("ci low") I`row'=("ci high")
	putexcel K`row'=("Model AIC")

	local row = 2
	local ++row
	putexcel A`row'=("Age <5 at SM episode")
	local ++row
	foreach var in mulout ecvt coatsc {
		qui run_mixed_analysis `var' "" "adjusted"
		write_results `row' `var'
		local ++row
	}
	local ++row
	putexcel A`row'=("Age <5 at SM episode; >=5 at follow-up cognitive testing")
	local ++row
	foreach var in mpi vdpt seqpro {
		qui run_mixed_analysis `var' "agecogt==2" "adjusted"
		write_results `row' `var'
		local ++row
	}
	local ++row
	putexcel A`row'=("Age >=5 at SM episode")
	local ++row
	foreach var in mpi vdpt seqpro {
		qui run_mixed_analysis `var' "agecogt==3" "adjusted"
		write_results `row' `var'
		local ++row
	}

	// unadjusted 
	putexcel set Fig5, modify sheet("5a_C1_Unadjusted")

	local row = 1
	putexcel C`row'=("N(obs), n") D`row'=("Beta, (95% CI)") E`row'=("P")
	putexcel F`row'=("Unadjusted") G`row'=("Beta") H`row'=("ci low") I`row'=("ci high")
	putexcel K`row'=("Model AIC")

	local row = 2
	local ++row
	putexcel A`row'=("Age <5 at SM episode")
	local ++row
	foreach var in mulout ecvt coatsc {
		qui run_mixed_analysis `var' "" "unadjusted"
		write_results `row' `var'
		local ++row
	}
	local ++row
	putexcel A`row'=("Age <5 at SM episode; >=5 at follow-up cognitive testing")
	local ++row
	foreach var in mpi vdpt seqpro {
		qui run_mixed_analysis `var' "agecogt==2" "unadjusted"
		write_results `row' `var'
		local ++row
	}
	local ++row
	putexcel A`row'=("Age >=5 at SM episode")
	local ++row
	foreach var in mpi vdpt seqpro {
		qui run_mixed_analysis `var' "agecogt==3" "unadjusted"
		write_results `row' `var'
		local ++row
	}

	
// Supp Table 5e

	preserve 
	keep if watson_sm==1
	
	// adjusted 
	putexcel set Supp_Tab5e, modify sheet("C1_Adjusted")

	local row = 1
	putexcel C`row'=("N(obs), n") D`row'=("Beta, (95% CI)") E`row'=("P")
	putexcel F`row'=("Adjusted") G`row'=("Beta") H`row'=("ci low") I`row'=("ci high")
	putexcel K`row'=("Model AIC")

	local row = 2
	local ++row
	putexcel A`row'=("Age <5 at SM episode")
	local ++row
	foreach var in mulout ecvt coatsc {
		qui run_mixed_analysis `var' "" "adjusted"
		write_results `row' `var'
		local ++row
	}
	local ++row
	putexcel A`row'=("Age <5 at SM episode; >=5 at follow-up cognitive testing")
	local ++row
	foreach var in mpi vdpt seqpro {
		qui run_mixed_analysis `var' "agecogt==2" "adjusted"
		write_results `row' `var'
		local ++row
	}
	local ++row
	putexcel A`row'=("Age >=5 at SM episode")
	local ++row
	foreach var in mpi vdpt seqpro {
		qui run_mixed_analysis `var' "agecogt==3" "adjusted"
		write_results `row' `var'
		local ++row
	}

	// unadjusted 
	putexcel set Supp_Tab5e, modify sheet("C1_Unadjusted")

	local row = 1
	putexcel C`row'=("N(obs), n") D`row'=("Beta, (95% CI)") E`row'=("P")
	putexcel F`row'=("Unadjusted") G`row'=("Beta") H`row'=("ci low") I`row'=("ci high")
	putexcel K`row'=("Model AIC")

	local row = 2
	local ++row
	putexcel A`row'=("Age <5 at SM episode")
	local ++row
	foreach var in mulout ecvt coatsc {
		qui run_mixed_analysis `var' "" "unadjusted"
		write_results `row' `var'
		local ++row
	}
	local ++row
	putexcel A`row'=("Age <5 at SM episode; >=5 at follow-up cognitive testing")
	local ++row
	foreach var in mpi vdpt seqpro {
		qui run_mixed_analysis `var' "agecogt==2" "unadjusted"
		write_results `row' `var'
		local ++row
	}
	local ++row
	putexcel A`row'=("Age >=5 at SM episode")
	local ++row
	foreach var in mpi vdpt seqpro {
		qui run_mixed_analysis `var' "agecogt==3" "unadjusted"
		write_results `row' `var'
		local ++row
	}

	restore

