setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(gtsummary)          
library(flextable)
library(openxlsx)
library(mice)
library(nhanesR)
library(pROC)
library(survival)
library(survminer)
source("3_zl_R_Function_v2.r")
load("IBD_data_ham_all.RData")
IBD_data_ham_all[,-1:-14] <- apply(IBD_data_ham_all[,-1:-14],2,as.numeric)

#mice多重插补
IBD_mice <- mice(IBD_data_ham_all,m=1,seed = 5)

#读取多重插补信息,并替换列名
IBD_mice_data <- complete(IBD_mice,action = 1)
IBD_mice_data$RAR <- round(IBD_mice_data$RDW/IBD_mice_data$ALB,2)
IBD_mice_data_nosur <- IBD_mice_data |> filter(IBD==TRUE) 
IBD_mice_data_nosur <- IBD_mice_data_nosur[,-11]

IBD_mice_data_nosur$svrvi_time <- ifelse(IBD_mice_data_nosur$svrvi_time=="no",0,
                                         IBD_mice_data_nosur$svrvi_time) |> as.numeric()

IBD_mice_data_nosur$RARQ <- quant(IBD_mice_data_nosur$RAR,round = 2,n=4,Q=T)
IBD_mice_data_dead <- IBD_mice_data_nosur |> filter(Statue==1)
IBD_mice_data_nosur$Inho_5 <- ifelse(IBD_mice_data_nosur$Inhospital>5,1,0)

#统计数据
tbl_summary(IBD_mice_data_nosur,by = "RARQ",
            statistic = list(all_continuous()~"({mean} ± {sd})",
                             all_categorical() ~ "{n}({p})"),
            digits = all_continuous() ~ 2) |> 
    add_p()|> 
    as_flex_table() |> 
    save_as_docx(path = "RARQ_IBD_diff.docx")

tbl_summary(IBD_mice_data_nosur,by = "Inho_5",
            statistic = list(all_continuous()~"({mean} ± {sd})",
                             all_categorical() ~ "{n}({p})"),
            digits = all_continuous() ~ 2) |> 
    add_p()|> 
    as_flex_table() |> 
    save_as_docx(path = "Inho_5_IBD_diff.docx")

#需要校正的数据
colnames(IBD_mice_data_nosur) |> dput()

adj_var <- c("RAR","ALB","RDW","RARQ")

adj_var_list <- name_list(adj_var)

#校正1，用年龄和性别
ad_mol1_index <- c("SEX","Age")

#校正2，用较正1+病史
ad_mol2_index <- c(ad_mol1_index,"nicotine","hyperlipidemia","hypertension","RA")

#校正3，用校正2+部分临床化验
ad_mol3_index <- c(ad_mol2_index, "Creatinine", 
 "WBC", "BUN", "Chloride", "Bicarbonate", "Phosphate", 
"Calcium", "INR", "PT", "AST", "ALT")

#adj模型0
adj_list <- list(NULL,ad_mol1_index,ad_mol2_index,ad_mol3_index)
lapply(adj_list,function(x) adj_log_glm(data=IBD_mice_data_nosur,
                   signle_list=adj_var_list,
                   adj_immo=x,
                   group_index = "Inho_5")) 

UC_data <- IBD_mice_data_nosur |> filter(UC==1)
lapply(adj_list,function(x) adj_log_glm(data=UC_data,
                                        signle_list=adj_var_list,
                                        adj_immo=x,
                                        group_index = "Inho_5"))

CD_data <- IBD_mice_data_nosur |> filter(CD==1) 
lapply(adj_list,function(x) adj_log_glm(data=CD_data,
                                        signle_list=adj_var_list,
                                        adj_immo=x,
                                        group_index = "Inho_5")) 

#分析随访人群
IBD_mice_data_nosur_svr <- IBD_mice_data_nosur |> filter(Statue==1)
fivenum(as.numeric(IBD_mice_data_nosur_svr$svrvi_time))
IBD_mice_data_nosur_svr$svrvi_time <- as.numeric(IBD_mice_data_nosur_svr$svrvi_time)

IBD_mice_data_nosur_svr$svr_group <- ifelse(IBD_mice_data_nosur_svr$svrvi_time>90,0,1)

tbl_summary(IBD_mice_data_nosur_svr,by = "svr_group",
            statistic = list(all_continuous()~"({mean} ± {sd})",
                             all_categorical() ~ "{n}({p})"),
            digits = all_continuous() ~ 2) |> 
    add_p() |> as_flex_table() |> 
    save_as_docx(path = "IBD_follow_diff.docx")

lapply(adj_list,function(x) adj_cox_reg(data=IBD_mice_data_nosur_svr,
            signle_list=adj_var_list,
            statue="svr_group",
            adj_immo=x,
            live_time="svrvi_time"))

#RCS寻找最佳的分界值
#寻找RAR的最佳RCS
ddist <- datadist(IBD_mice_data_nosur_svr)
options(datadist="ddist")
serch_index <- "ALB"
adju_index <-  "Age"

S <- Surv(IBD_mice_data_nosur_svr$svrvi_time,IBD_mice_data_nosur_svr$svr_group==1)
for (knot in 3:10) {
    cph_formula <- as.formula(paste0("S~rcs(",serch_index,",",knot,")+",adju_index))
    
    fit <- cph(cph_formula,data=IBD_mice_data_nosur_svr,x= TRUE, y= TRUE, surv = TRUE)
    tmp <- extractAIC(fit)
    if(knot==3){AIC=tmp[2];nk=3}
    if(tmp[2]<AIC){AIC=tmp[2];nk=knot}
}


#---------基于cox绘制 单因素LnAl 的 RCS
pacman::p_load(rms,survminer,ggplot2,ggsci)
#构建cph 函数获取rcs, 单指标LnAl；也可校正其他因素 S ~ rcs(LnAl,4)+ BMI 等
cph_formula <- as.formula(paste0("S~rcs(",serch_index,",",nk,")+",adju_index))

fit.RAR <- cph(cph_formula, x=TRUE, y=TRUE,data=IBD_mice_data_nosur_svr)
# PH 检验,P>0.05 符合PH假设
cox.zph(fit.RAR, "rank")    
#残差图的横轴是时间，纵轴是残差，残差均匀分布则表示残差与时间相互独立，满足ph假设
ggcoxzph(cox.zph(fit.RAR, "rank"))

pdf_name <- paste0("RCS_",serch_index,"-",adju_index,".pdf")
pdf(pdf_name,width = 5,height = 5)
ggcoxzph(cox.zph(fit.RAR, "rank")) 
dev.off()

#非线性检验p-non-linear，0.9289
anova(fit.RAR)                      
p <-round(anova(fit.RAR)[,3],3)
## HR计算，fun是转化函数
Pre_HR.RAR <-rms::Predict(fit.RAR,ALB,fun=exp,type="predictions",ref.zero=T,conf.int = 0.95,digits=2)
ggplot(Pre_HR.RAR)
# y-hat=1.00000, 对应LnAl= 5.006
#View(Pre_HR.RAR)

# 全人群的COX模型 RCS绘图
rcs_all <-  ggplot()+
    geom_line(data=Pre_HR.RAR,
              aes(ALB,yhat),# x y轴距
              linetype="solid",#曲线加粗
              size=1,
              alpha=0.7,
              colour="green")+
    scale_fill_nejm()+ ##采用ggsci包中英格兰调色，也可以其他
    geom_ribbon(data=Pre_HR.RAR,
                aes(ALB, ymin=lower,ymax=upper,fill="zx"),alpha=0.1)+
    theme_classic()+
    scale_fill_nejm()+
    geom_hline(yintercept=1,linetype=2,size=0.75) #y=1水平线+
labs(title ="风险随LnAl变化曲RCS",
     x="ALB", 
     y="HR (95%CI)"
)
pdf_name1 <- paste0("rcs_",serch_index,"_",adju_index,"_.pdf")
pdf(pdf_name1,width = 5,height = 4)
rcs_all
dev.off()

#KM曲线
IBD_mice_data_nosur_svr$RAR_5.64 <- ifelse(IBD_mice_data_nosur_svr$RAR>5.64,"high","low")
IBD_mice_data_nosur_svr$RDW_15.9 <- ifelse(IBD_mice_data_nosur_svr$RDW>15.9,"high","low")
IBD_mice_data_nosur_svr$ALB_2.9 <- ifelse(IBD_mice_data_nosur_svr$ALB>2.9,"high","low")

#KM分析的批量做图
IBD_mice_data_nosur_svr1 <- IBD_mice_data_nosur_svr
IBD_mice_data_nosur_svr1$svrvi_time <- ifelse(IBD_mice_data_nosur_svr1$svrvi_time>120,120,IBD_mice_data_nosur_svr1$svrvi_time)


Var_KM <- grep("_KM",colnames(IBD_mice_data_nosur_svr1),value = T)
Var_KM <- c("RAR_5.64","RDW_15.9","ALB_2.9")
km_formula <- sapply(Var_KM,function(x) as.formula(paste0('Surv(svrvi_time, svr_group)~',x)))
km_analys <- lapply(km_formula,function(x) survfit(x,data=IBD_mice_data_nosur_svr1))

for(i in Var_KM){
    
    fit_formula <- as.formula(paste0('Surv(svrvi_time, svr_group)~',i))
    fit <- survfit(fit_formula,data=IBD_mice_data_nosur_svr1)
    fit$call$formula <- fit_formula
    
    fig <- ggsurvplot(fit,
                      pval = TRUE, 
                      pval.coord=c(0.5,0.7),
                      conf.int = TRUE,
                      risk.table = TRUE, 
                      risk.table.col = "strata", 
                      linetype = "strata", #fun = "cumhaz",
                      surv.median.line = "hv", # 同时显示垂直和水平参考线
                      ggtheme = theme_bw(), 
                      xlab = "Follow up time(d)",
                      #legend = c(0.8,0.75), 
                      legend.title = "",
                      legend.labs = c("High", "Low"),
                      break.x.by = 20,
                      ylim = c(0.5,1),
                      palette = "hue")
    
    pdf(file = paste0(i,"_mimic.pdf"),width=8, height=6)
    print(fig,newpage = FALSE)
    dev.off()
}



#构建列线图
#将数据打包好
IBD_mice_data_nosur_svr$RARQ <- quant(IBD_mice_data_nosur_svr$RAR, n = 4,Q = TRUE,round=5)
IBD_mice_data_nosur_svr$RARQ.median <- quant.median(IBD_mice_data_nosur_svr$RAR, n = 4,round=2)

ddist <- datadist(IBD_mice_data_nosur_svr)
options(datadist='ddist')

#构建多因素的Cox回归模型
cox <- cph(Surv(svrvi_time, svr_group) ~ Age+RARQ.median+ALB+RDW,
           data = IBD_mice_data_nosur_svr,x=T,y=T,surv = T)

surv <- Survival(cox)
sur_1_year<-function(x)surv(365,lp=x)#1年生存
sur_2_year<-function(x)surv(730,lp=x)#2年生存
sur_3_year<-function(x)surv(1095,lp=x)#3年生存
#sur_4_year<-function(x)surv(1460,lp=x)#4年生存
#sur_5_year<-function(x)surv(1825,lp=x)#5年生存

# 做列线图
nom_sur <- nomogram(cox,fun=list(sur_1_year),
                    lp= F,
                    funlabel=c('1-Year Survival'),maxscale=100,fun.at=c('0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1'))

pdf("nomogram.pdf",width=8, height=6)
plot(nom_sur)
dev.off()