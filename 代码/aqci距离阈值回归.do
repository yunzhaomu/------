clear all
global y "z_aqci"
global process_x "aqci lndist ET rain temperature elevation slope light pop road_density"
global thresh_x  "z_ET z_rain z_temperature z_elevation z_slope z_light z_pop z_road_density"
global vars "r1_lndist r2_lndist ET rain temp elev slope light pop road"

import delimited "D:\毕业论文数据\导出表格数据\RS_park_mountain.csv", case(preserve) clear

* --- 计算空气质量分指标并合成总指标
foreach var in AOD SO2 CO O3 NO2 {
    * 计算分位数并存入标量 (逻辑保持不变)
    centile `var', centile(20 40 60 80)
    scalar b1_`var' = r(c_1)
    scalar b2_`var' = r(c_2)
    scalar b3_`var' = r(c_3)
    scalar b4_`var' = r(c_4)
    
    summarize `var'
    scalar b0_`var' = r(min)
    scalar b5_`var' = r(max)
}

foreach var in AOD SO2 CO O3 NO2 {
    gen iaqi_`var' = .
    replace iaqi_`var' = (50 - 0) / (b1_`var' - b0_`var') * (`var' - b0_`var') + 0 ///
        if `var' >= b0_`var' & `var' <= b1_`var'
        
    replace iaqi_`var' = (100 - 50) / (b2_`var' - b1_`var') * (`var' - b1_`var') + 50 ///
        if `var' > b1_`var' & `var' <= b2_`var'
        
    replace iaqi_`var' = (150 - 100) / (b3_`var' - b2_`var') * (`var' - b2_`var') + 100 ///
        if `var' > b2_`var' & `var' <= b3_`var'
        
    replace iaqi_`var' = (200 - 150) / (b4_`var' - b3_`var') * (`var' - b3_`var') + 150 ///
        if `var' > b3_`var' & `var' <= b4_`var'
        
    replace iaqi_`var' = (300 - 200) / (b5_`var' - b4_`var') * (`var' - b4_`var') + 200 ///
        if `var' > b4_`var'
}
egen aqci = rowmean(iaqi_AOD iaqi_SO2 iaqi_CO iaqi_O3 iaqi_NO2)


levelsof year, local(years)
gen keep_mark = 0
egen id = group(coordinate)
duplicates drop id year,force
gen lndist = ln(dist_ring + 1)
foreach v in $process_x {
    egen z_`v' = std(`v')
}

* 筛选控制变量集
// 检验多重共线性
reg $y lndist $thresh_x
vif
// 逐步回归筛选
sw reg $y lndist $thresh_x, pe(.05)

* --- 准备工作 ---
tempname memhold
* 定义存储文件：包括门限、两个区间的核心变量、以及所有控制变量
postfile `memhold' obs thr $vars using "mountain_aqci_results.dta", replace

* --- 循环迭代 ---
forvalues i = 1/20 {
	foreach j in `years' {
		preserve
		keep if year == `j'
		count if type == "inner"
		local inner_cnt = r(N)

		count if type == "buffer"
		local buffer_cnt = r(N)
		local sample_cnt : display ceil(min(2500,`inner_cnt')*3/40)

		set seed `i'  // 设置种子保证结果可重复
		gen rand_num = runiform() // 生成随机数

		* 1. 保留2000个内部区点 (dist_ring == 0)和大概3倍的内部区点
		bysort dist_ring (rand_num): gen ring_count = _n
		replace keep_mark = 1 if type == "inner" & ring_count <= 2500
		replace keep_mark = 1 if type == "buffer" & ring_count <= `sample_cnt'

		* 3. 只保留选中的点进行后续分析
		keep if keep_mark == 1
		
		* 运行门限模型
		quietly threshold $y $thresh_x, ///
			threshvar(dist_ring) regionvars(lndist) nthresholds(1) vce(robust)
		
		* 样本数
		local n_val = _N
		* 提取门限值
		matrix mat_thr = e(thresholds)
		local t_val = mat_thr[1,2]
		* 提取核心变量系数 (两个区间)
		local b1 = _b[Region1:lndist]
		local b2 = _b[Region2:lndist]
		if _rc == 0 & `t_val' > 1{			
			* 提取区间不变的控制变量系数
			local b_et   = _b[z_ET]
			local b_rain = _b[z_rain]
			local b_temp = _b[z_temperature]
			local b_elev = _b[z_elevation]
			local b_slop = _b[z_slope]
			local b_ligh = _b[z_light]
			local b_pop = _b[z_pop]
			local b_road = _b[z_road_density]
			
			post `memhold' (`n_val') (`t_val') (`b1') (`b2') (`b_et') (`b_rain') (`b_temp') (`b_elev') (`b_slop') (`b_ligh') (`b_pop') (`b_road')
		}
		restore
    }
    if mod(`i', 10) == 0 noisily display "Iteration `i' finished"
}
postclose `memhold'

* --- 读取 1000 次模拟的结果文件 ---
use "mountain_aqci_results.dta", clear
local vars "r1_lndist r2_lndist ET rain temp elev slope light pop road"

* --- 计算均值 (系数) 和 标准差 (标准误) ---
* 提取均值到矩阵 b
tabstat `vars', statistics(mean) save
matrix b = r(StatTotal)
local b_r1 = b[1, "r1_lndist"]
local b_r2 = b[1, "r2_lndist"]

* 提取标准差到矩阵 SE，并转化为方差矩阵 V (对角阵)
tabstat `vars', statistics(sd) save
matrix se_diag = r(StatTotal)
matrix V = diag(se_diag) // 将标准误转为方差 (方差=标准误^2)，但outreg2通常通过矩阵对角线识别
* 注意：V矩阵的维度必须是 k*k，这里构造一个对角阵
forvalues j = 1/`: word count `vars'' {
    matrix V[`j',`j'] = se_diag[1,`j']^2
}

* --- 3. 设置矩阵的行列名称 (必须与后续变量名匹配) ---
matrix colnames b = `vars'
matrix colnames V = `vars'
matrix rownames V = `vars'

* 计算平均门限值
quietly sum thr
local mean_thr : display %9.1f r(mean)
quietly sum obs
local obs : display r(mean)

* --- 4. 伪造回归环境 ---
ereturn post b V,obs(`obs') depname($y)

* --- 5. 使用 outreg2 导出到 Word/Excel ---
outreg2 using "D:\毕业论文数据\论文图表\mountain_aqci_Bootstrap_Result.doc", replace ///
    word bdec(3) tdec(3) rdec(3) ///
    title("Table: Threshold Model with 1000 Bootstrap Replications") ///
    addnote("Threshold variable: dist_ring. Controlled variables are standardized.") ///
    addtext("Threshold","`mean_thr'","Bootstrap Reps", "1000")

* 2. 打开 Excel 准备写入
putexcel set "D:\毕业论文数据\结果\mountain_aqci_Bootstrap_Result.xlsx", replace

* 3. 写入表头和数值 (A列为名称，B列为数值)
putexcel A1 = ("Parameter") B1 = ("Mean Value")
putexcel A2 = ("Threshold")  B2 = (`mean_thr')
putexcel A3 = ("r1_lndist") B3 = (`b_r1')
putexcel A4 = ("r2_lndist") B4 = (`b_r2')

* 4. 格式化（可选：设置加粗和数字格式）
putexcel A1:B1, bold border(bottom)
putexcel B2, nformat("0")
putexcel B2:B4, nformat("0.000")

// xtset id year
// xthreg z_aqci z_ET z_rain z_temperature z_light, rx(lndist) qx(dist_ring) thnum(1) fe trim(0.5) bs(1000) // 这个要求time-varying