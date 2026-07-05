clear all
set seed 24
global park "rainforest"
global y "CO"
if "$park" == "pandas" {
    global thr 12.6
}
else if "$park" == "tiger" {
    global thr 15.8
}
else if "$park" == "rivers" {
    global thr 3.9
}
else if "$park" == "rainforest" {
    global thr 6.8
}
else if "$park" == "mountain" {
    global thr 12
}
global psm_control  "rain temperature elevation light"

import delimited "D:\毕业论文数据\导出表格数据\RS_park_${park}.csv", case(preserve) clear

* 1. 环境准备与数据清洗
gen treat = (dist_ring >0 & dist_ring <= $thr)
drop if dist_ring == 0 // 排除内部区

* 2. 定义存放结果的 postfile
* 增加 red_mbias (平均偏差减少率) 和 red_medbias (中位数偏差减少率)
tempname results
postfile `results' year mean_treat mean_ctrl att pseudoR2_pre pseudoR2_post lr_pre lr_post ///
                  m_bias_pre m_bias_post red_mbias med_bias_pre med_bias_post red_medbias b_pre b_post ///
                  using "${park}_${y}_buffer_control.dta", replace

* 3. 年份循环
levelsof year, local(years)

foreach y in `years' {
    preserve
    keep if year == `y'
    
    * --- 匹配执行 ---
	psmatch2 treat $psm_control , outcome($y) logit neighbor(5) ties common ate caliper(0.05)
    
    * 提取均值和 ATT
    summarize $y if treat == 1 & _support == 1
    local m1 = r(mean)
    summarize $y [aw=_weight] if treat == 0 & _support == 1
    local m0 = r(mean)
    local att = `m1' - `m0'
    
    * --- 提取平衡性指标 ---
    * 使用 pstest 获取匹配前后的全局统计量
    pstest $psm_control, both
    
    * Pseudo R2 & LR Chi2
    local r2_pre   = r(r2bef)
    local r2_post  = r(r2aft)
    local lr_pre    = r(chiprobbef)
    local lr_post   = r(chiprobaft)
    
    * Mean Bias & 中位数 Bias
    local mb_pre    = r(meanbiasbef)
    local mb_post   = r(meanbiasaft)
    local medb_pre  = r(medbiasbef)
    local medb_post = r(medbiasaft)
    
    * 计算减少的百分比 (Reduction %)
    local r_mbias   = (`mb_pre' - `mb_post') / `mb_pre' * 100
    local r_medbias = (`medb_pre' - `medb_post') / `medb_pre' * 100
    
    * Rubin's B
    local b_pre     = r(Bbef)
	local b_post     = r(Baft)
    
    * --- 写入数据 ---
    post `results' (`y') (`m1') (`m0') (`att') (`r2_pre') (`r2_post') (`lr_pre') (`lr_post') ///
                  (`mb_pre') (`mb_post') (`r_mbias') (`medb_pre') (`medb_post') (`r_medbias') (`b_pre') (`b_post')
    
    di "年份 `y' 统计完成。"
    restore
}

postclose `results'

* 4. 载入结果并格式化展示
use "${park}_${y}_buffer_control.dta", clear

* 最终列表预览
list year att m_bias_pre m_bias_post red_mbias med_bias_pre med_bias_post red_medbias b_pre b_post, clean

export excel using "D:/毕业论文数据/结果/${park}_${y}_buffer_control.xlsx", firstrow(variables) replace
