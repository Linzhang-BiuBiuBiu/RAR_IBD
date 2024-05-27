library(dplyr)
library(fBasics)
library(readxl)

#计算NA值的函数，data是数据集，number为1是行选择，2为列
na_count <- function(data,number){
    na_count_result <- apply(data,number,function(x) sum(is.na(x)))
    print(na_count_result)
}

#删除NA值过多的列，data为数据集，number为NA的阈值，超过阈值则删除
na_delet_data <- function(data,number){
    for(i in colnames(data)){
        na_count_tem <- sum(is.na(data[,i]))
        if(na_count_tem >= number){
            data <- data[,which(colnames(data)!= i)]
        }
    }
    return(data)
}

#删除向量函数
drop_vector <- function(x,y){
    for(i in y){
        if(i %in% x){
            where_i <- which(x %in% i)
            x <- x[-where_i]
        }
    }
    return(x)
}

#缺失值用中位数进行填补
na_median_sub <- function(analysis_data,number){
for(i in colnames(analysis_data)){
    na_count <- sum(is.na(analysis_data[,i]))
    if (na_count > number) {
        analysis_data <- analysis_data[,colnames(analysis_data)!=i]
    }else{
        na_mean <- median(as.numeric(analysis_data[!is.na(analysis_data[,i]),i]))
        where_na <- is.na(analysis_data[,i])
        analysis_data[where_na,i] <- na_mean}}
   return(analysis_data)}

#去掉>,<,--,中文字符，非数字符号
del_sign <- function(data){
    for(i in colnames(data)){
        data[,i] <- gsub("<|>|-","",data[,i])
    }
    return(data)
}

#差异分析,默认group在第一列
differenc_test <- function(data,group){
    colname_data <- drop_vector(colnames(data),group)
    chi_t_data <- c()
    for (i in colname_data) {
        unique_numb <- na.omit(unique(data[,i]))
        if(length(unique_numb) < 3){
        data_i <- data[,c(i,group)]
        data_chitest <- data_i[!is.na(data_i[,i]),]
        chitest <- chisq.test(data_chitest[,1],data_chitest[,2])
        chi_p <- round(chitest$p.value,4)
        sum_chi0 <- table(data_chitest[data_chitest[,2]==0,i])
        sum_chi1 <- table(data_chitest[data_chitest[,2]==1,i])
        
        sum_pre0 <- round(sum_chi0[2]/sum(sum_chi0),4)*100
        sum_pre1 <- round(sum_chi1[2]/sum(sum_chi1),4)*100
        
        chi_all0 <- paste0(sum_chi0[2],"(",sum_pre0,")")
        chi_all1 <- paste0(sum_chi1[2],"(",sum_pre1,")")
        
        diff_matrix <- c(value=i,distrbution0=chi_all0,distrbution1=chi_all1,P_value=chi_p)
        
    }else{
        tformula<- as.formula(paste0(i,"~",group))
        data_na_no <- data[!is.na(data[,i]),c(i,group)]
        ttest <- t.test(tformula,data_na_no)
        test_p <- round(ttest$p.value,4)
        mean0 <- round(ttest$estimate[1],2)
        mean1 <- round(ttest$estimate[2],2)
        
        sd0 <- round(sd(data_na_no[data_na_no[,2]==0,1]),2)
        sd1 <- round(sd(data_na_no[data_na_no[,2]==1,1]),2)
        
        mean_sd0 <- paste0(mean0," ± ",sd0)
        mean_sd1 <- paste0(mean1," ± ",sd1)
        
        diff_matrix <- c(value=i,distrbution0=mean_sd0,distrbution1=mean_sd1,P_value=test_p)
    }
    chi_t_data <- rbind(diff_matrix,chi_t_data)}   
    return(chi_t_data)
}


#批量单因素logistic
uni_log_test <- function(data,group){
    data <- as.data.frame(data)
    colname_data <- drop_vector(colnames(data),group)
    uni_glm_result <- c()
    for(i in colname_data){
        formula_glm <- formula(paste0(group,"~",i))
        glm1<-glm(formula = formula_glm,data=data,family= "binomial")  
        
        glm2<-summary(glm1)
        
        OR<-round(exp(coef(glm1)[-1]),3)
        
        SE<-glm2$coefficients[,2][-1]
        
        CI5<-round(exp(coef(glm1)[-1]-1.96*SE),3) #95%CI
        
        CI95<-round(exp(coef(glm1)[-1]+1.96*SE),3)
        
        CI<-paste0(CI5,"-",CI95)
        
        P<-glm2$coefficients[,4][-1]
        
        uni_glm_model <- data.frame("Charactersitc"=rownames(as.data.frame(OR)),
                                    "OR"=OR,
                                    "CI"=CI,
                                    "P"=P)
        uni_glm_result <- rbind(uni_glm_result,uni_glm_model)
    }
    return(uni_glm_result)
}

#校正logistic
adj_log_test <- function(data,group,del_vec){
    data <- as.data.frame(data)
    colname_data <- drop_vector(colnames(data),group)
    colname_data <- drop_vector(colname_data,del_vec)
    
    foumula_glm <- formula(paste0(group,"~."))
    glm1 <- glm(foumula_glm,data[,c(colname_data,group)],family=binomial)
    glm2<-summary(glm1)
    
    OR<-round(exp(coef(glm1)[-1]),3)
    SE<-glm2$coefficients[,2][-1]
    
    CI5<-round(exp(coef(glm1)[-1]-1.96*SE),3) #95%CI
    CI95<-round(exp(coef(glm1)[-1]+1.96*SE),3)
    CI<-paste0(CI5,"-",CI95)
    
    P<-round(glm2$coefficients[,4][-1],5)
    
    mul_glm_model_RA <- data.frame("OR"=OR,
                                   "CI"=CI,
                                   "P"=P)
}

#单因素cox
uni_cox_test <- function(data,statue,followtime){
    con_var <- drop_vector(colnames(data),c(statue,followtime))
    tdsinglcox <- c()
    for(i in con_var){
        uni_cox_formula  <- as.formula(paste0("Surv(","liv_time",",","EXPIRE_FLAG",")~",i))
        
        coxpha <- coxph(uni_cox_formula, data = data)
        
        fit <- summary(coxpha)
        
        P<-signif(fit$coefficients[,5], digits=3)
        HR <-signif(fit$coef[,2], digits=3)
        
        HR.confint.lower <- signif(fit$conf.int[,"lower .95"], 3)
        HR.confint.upper <- signif(fit$conf.int[,"upper .95"],3)
        HR <- paste0(HR, " (",HR.confint.lower, "-", HR.confint.upper, ")")
        
        if(length(HR.confint.lower)>1){
            sig_cox <-data.frame("chara"=rownames(as.data.frame(HR.confint.lower)),"HR (95% CI for HR)"=HR,"p"=P)
        }else{
            sig_cox <-data.frame("chara"=i,"HR (95% CI for HR)"=HR,"p"=P)
        }
        tdsinglcox <- rbind(sig_cox,tdsinglcox)
    }
    return(tdsinglcox)
    
}


#多因素cox
adj_cox_test <- function(data,statue,followtime,del_var){
    con_var <- drop_vector(colnames(data),c(statue,followtime,del_var))
    
    adj_cox_formula  <- as.formula(paste0("Surv(",followtime,",",statue,")~."))
    tdmultiCox=coxph(adj_cox_formula,data[,c(con_var,statue,followtime)])
    
    tdmultiCoxSum=summary(tdmultiCox)
    
    outResult=cbind(
        HR=round(tdmultiCoxSum$conf.int[,"exp(coef)"],3),
        "HR (95% CI for HR)"=paste0(round(tdmultiCoxSum$conf.int[,"lower .95"],3),"-",round(tdmultiCoxSum$conf.int[,"upper .95"],3)),
        pvalue=round(tdmultiCoxSum$coefficients[,"Pr(>|z|)"],3))
    
    adj_cox_Result=cbind(id=row.names(outResult),outResult)
    
    return(adj_cox_Result)
}

#分段死亡率
quantile_mortality <- function(data,statue,varia,number){
    data <- data[!is.na(data[,varia]),]
    data <- as.data.frame(data)
    
    data[,"variabQ"] <- quant(data[,varia], n = number,Q = TRUE,round=4)
    data[,"variabQ.median"] <- quant.median(data[,varia], n = number,round=4)
    
    deth_rar_q <- c()
    for (i in unique(data[,"variabQ"])) {
        mordi_data <-data[data[,"variabQ"]==i,]
        dedth_indi <- mordi_data[mordi_data[,statue]==1,]
        deth_perc <- nrow(dedth_indi)/nrow(mordi_data)
        var_mor <- paste0(varia,"_mor")
        q_deth <- c(Q=i,var_mor=deth_perc)
        deth_rar_q <- rbind(q_deth,deth_rar_q)
    }
    return(deth_rar_q)
}

#为了获取患者的所有化验，把门诊号统一成住院号
#也就是把220510131778替换成00586417
#data是需要替换的数据，index_col是需要替换的列，index_data是构建索引的数据框
#primary是index_data需要被替换的行，current是index_data中替换后的行

sub_data_stri <- function(data=data,
                          index_col=index_col,
                          index_data=index_data,
                          primary=primary,
                          current=current){
    data <- as.data.frame(data)
    index_data_no_na <- na.omit(index_data[,c(primary,current)])
    for(i in 1:nrow(index_data_no_na)){
        
        data[data[,index_col] %in% index_data_no_na[i,primary],index_col] <- index_data_no_na[i,current] 
        
    }
    
    return(data)
}

#读取目前文件夹下的所有xlsx中的第x个sheet，并自动merge
merge_xlsx <- function(x){
    xlsx_file_names <- grep("xlsx$", dir(), value=T)
    xlsx_all <- data.frame()
    
    for(i in xlsx_file_names){
        tem_xlsx <- read_xlsx(i,x) 
        xlsx_all <- rbind(xlsx_all,tem_xlsx)
    }
    return(xlsx_all)
}

#读取列表中的第n个文件，并整合成data.frame
summary_list_n <- function(n,list=list_file){
    list_vector <- data.frame() 
    for (i in 1:length(list_file)) {
        tem_list <- list_file[[i]]
        list_vector <- rbind(list_vector,tem_list[n])
    }
    return(list_vector)
}

#提取患者的化验，data是化验集合，
#human_ID是住院号类的个人ID
#lab_col是化验名称
#lab_result是化验结果
#提取第n次化验
extrade_n_lab <- function(data=data,n=n,
                              human_ID = human_ID,
                              lab_col=lab_col,
                              lab_result=lab_result){
    
    data <- as.data.frame(data)
    
    data <- data[!is.na(data[,lab_col]),]
    data <- data[!is.na(data[,lab_result]),]
    data <- data[!is.na(data[,human_ID]),]
    
    lab_names <- t(as.data.frame(unique(data[,lab_col])))
    ID_lab_names <- cbind(human_ID,lab_names)
    colnames(ID_lab_names) <- ID_lab_names[1,]
    
    unique_human_ID <- na.omit(as.data.frame(unique(data[,human_ID])))
    
    for(i in unique_human_ID[,1]){
      
      print(which(unique_human_ID[,1] %in% i))
        tem_uni_lab <- unique(data[data[,human_ID] %in% i,])
        
        tem_lab_result <- c()
        for(j in colnames(ID_lab_names)[-1]){
            
            try(lab_tem_resu <- tem_uni_lab[tem_uni_lab[,lab_col] %in% j,][n,],T)
            
            if(nrow(lab_tem_resu)>0){
                tem_lab_result <- rbind(tem_lab_result,lab_tem_resu[,c(lab_col,lab_result)])
            }
        }
        tem_lab_result <- na.omit(tem_lab_result)
        if(nrow(tem_lab_result)>0){
        tem_lab_result <- rbind(human_ID,tem_lab_result)
        rownames(tem_lab_result) <- tem_lab_result[,1]
        colnames(tem_lab_result)[2] <- i
        tem_lab_result[1,2] <- i
        tem_lab_result <- t(tem_lab_result)
        tem_lab_result <- t(tem_lab_result[-1,])
        ID_lab_names <- merge(tem_lab_result,ID_lab_names,all=T)
        if(ID_lab_names[2,1]==human_ID){
            ID_lab_names <- ID_lab_names[-2,] 
        }}}
    return(ID_lab_names)
} 

#提取最大的化验
extrade_max_lab <- function(data=data,
                              human_ID = human_ID,
                              lab_col=lab_col,
                              lab_result=lab_result){
  
  data <- as.data.frame(data)
  
  data <- data[!is.na(data[,lab_col]),]
  data <- data[!is.na(data[,lab_result]),]
  data <- data[!is.na(data[,human_ID]),]
  
  data[is.na(data)]<-0
  lab_names <- t(as.data.frame(unique(data[,lab_col])))
  ID_lab_names <- cbind(human_ID,lab_names)
  colnames(ID_lab_names) <- ID_lab_names[1,]
  
  unique_human_ID <- na.omit(as.data.frame(unique(data[,human_ID])))
  
  for(i in unique_human_ID[,1]){
    tem_uni_lab <- unique(data[data[,human_ID] %in% i,])
    
    tem_lab_result <- c()
    for(j in colnames(ID_lab_names)[-1]){
      
      try(lab_tem <- tem_uni_lab[tem_uni_lab[,lab_col] %in% j,],T)
      
      if(nrow(lab_tem)>0){
        lab_tem_resu <- lab_tem[which(lab_tem[,lab_result] %in% max(lab_tem[,lab_result])),][1,]
      }else{
        lab_tem_resu <- lab_tem[1,]
      }
      
      if(nrow(lab_tem_resu)>0){
        tem_lab_result <- rbind(tem_lab_result,lab_tem_resu[,c(lab_col,lab_result)])
      }
    }
    tem_lab_result <- na.omit(tem_lab_result)
    tem_lab_result <- rbind(human_ID,tem_lab_result)
    rownames(tem_lab_result) <- tem_lab_result[,1]
    colnames(tem_lab_result)[2] <- i
    tem_lab_result[1,2] <- i
    tem_lab_result <- t(tem_lab_result)
    tem_lab_result <- t(tem_lab_result[-1,])
    ID_lab_names <- merge(tem_lab_result,ID_lab_names,all=T)
    if(ID_lab_names[2,1]==human_ID){
      ID_lab_names <- ID_lab_names[-2,] 
    }}
  return(ID_lab_names)
} 

#读取对应目录下的所有excel的第n个sheet，并做成list
xlsx_list <- function(wd=wd,n=1){
    xlsx_files <-  grep("xlsx$", dir(wd), value=T)
    
    xlsx_list_file <- list()
    for (i in xlsx_files) {
        try(tem_xlsx <- read_xlsx(i,n),T)
        if(!is.null(tem_xlsx)&nrow(tem_xlsx)>0){
            xlsx_list_file[[i]] <- tem_xlsx
        }
    }
    return(xlsx_list_file)
}

#从左到右，从上到下，将字符串串联（用来提取手术记录的结果），根据首尾特征截取字符串，
factor_combine <- function(factor_data=factor_data,start1="结论",end1="术后医嘱"){
    row_data_frame <- data.frame()
    for(i in 1:nrow(factor_data)){
       tem_row_data <- paste(factor_data[i,],collapse = "")
       tem_row_data <- gsub("NA","",tem_row_data)
       row_data_frame <- rbind(row_data_frame,tem_row_data)
    }
    colnames(row_data_frame) <- "tem"
    factor_combine_data <- paste(row_data_frame,collapse = "")
    factor_cut1 <- strsplit(factor_combine_data,split=start1)[[1]][2]
    factor_cut2 <- strsplit(factor_cut1,split=end1)[[1]][1]
    return(factor_cut2)
}

#脂质组的PCA函数，其中列是脂质名称，file是输出的pdf名字，main是标题名字
PCA_lipid <- function(data=data,
                      nor=1,
                      file=file,
                      main=main){
    data[is.na(data)]<- 0
    data <- data[rowMeans(data)>0,]
    data <- data[,colMeans(data)>0]
    
    rowna_data <- rownames(data)
    if(nor==1){
        data <- apply(data, 1, function(x) x/sum(x))
        data.pca=prcomp(t(data))
    }else{data.pca=prcomp(data)}
    
    
    pcaPredict=predict(data.pca)
    
    type_data <- strsplit2(rowna_data,"-")[,2]
    
    PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=type_data)
    
    pdf(file=file, width=5.5, height=4.25)
    p1=ggscatter(data=PCA, x="PC1", y="PC2", color="Type", shape=19, 
                 ellipse=F, ellipse.type="norm", ellipse.border.remove=F, ellipse.alpha = 0.1,
                 size=2, main=main, legend="right")+
        theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))
    print(p1)
    dev.off()
    
}

#QC校正方法，使用QC来校正数据，目前只能用adj==1
qc_normalize <- function(primar_data=primar_data,adj_par=1){
    #定位QC在第几针
    
    Type <- strsplit2(colnames(primar_data),"-")[,3]
    
    #如果该序列的最后一针是QC，那就用后面的QC来校正前面的数据
    qc_index <- which(Type=="QC")
    
    mean_qc <- as.data.frame(apply(primar_data[,qc_index],1,mean))
    colnames(mean_qc) <- paste0(paste(strsplit2(colnames(primar_data),"-")[1,c(1,2)],collapse = "-"),"-QC-",max(qc_index)*sample(4:500,1))
        
    #如果最后没有以QC结尾，主动添加名称为QC—mean的列
    if(Type[ncol(primar_data)]!="QC"){
        primar_data1 <- cbind(primar_data,mean_qc)
    }else{primar_data1 <- primar_data}

    #重新注释type
    Type <- strsplit2(colnames(primar_data1),"-")[,3]
    
    #用前一针的QC校正后面的数据
    qc_index <- which(Type=="QC")
    
    qc_index1 <- qc_index[-1]
    qc_index2 <- qc_index[-length(qc_index)]
    
    qc_minus <- qc_index1-qc_index2
    qc_no_con_from <- qc_index[qc_minus!=1]
    qc_no_con_conut <- qc_minus[qc_minus!=1]-1
    
    qc_no_con_1 <- qc_index[qc_minus==1]
    #校正方法1，用QC的值除以每一针样品的值
    if(adj_par==1){
        for(i in 1:length(qc_no_con_from)){
            if(i ==1){
                adj_col_index <- 1:(qc_no_con_from[i]-1)
                primar_data1[,adj_col_index] <- primar_data1[,adj_col_index]/primar_data1[,(qc_no_con_from[i]-1)]   
            }
            adj_col_index <- qc_no_con_from[i]:(qc_no_con_from[i]+qc_no_con_conut[i])
            primar_data1[,adj_col_index] <- primar_data1[,adj_col_index]/primar_data1[,qc_no_con_from[i]]
        }
        primar_data1[,ncol(primar_data1)] <- primar_data1[,ncol(primar_data1)]/primar_data1[,ncol(primar_data1)]
        primar_data1[primar_data1=="NaN"] <- 0
        primar_data1[primar_data1=="Inf"] <- 0
        
        
    }

    #校正方法2，假设每一针QC是线性且样品间的数据梯度变化，以线性函数寻找值。
    if(adj_par==2){
        for(ind in 1:(length(qc_no_con_from)-1)){
            i <- qc_no_con_from[ind]
            j <- qc_no_con_from[ind+1]
            formula_line <- (primar_data1[j]/primar_data1[i])/(j-i)
            
            ck <- 1
            step_for <- formula_line
            while(ck+1 < (j-i)){
                ck <- ck+1
                step_for <- cbind(step_for,ck*formula_line)
            }
            primar_data1[,(i+1):(j-1)] <- primar_data1[,(i+1):(j-1)]/step_for
        }
    }
    return(primar_data1)
}

#合并脂质组得正负离子数据函数，输入形式首先是list，每个list得格式是一样的
lipid_data_merge <- function(lipid_list=lipid_list){
for(i in 1:length(names(lipid_list))){
    
    name_index <- names(lipid_list)[i]
    
    tem_file <- t(lipid_list[[name_index]])
    tem_file <- cbind(ID=rownames(tem_file),tem_file)
    rownames(tem_file) <- NULL
    if (i == 1) {
        data_all <- tem_file
    }else{data_all <- merge(data_all,tem_file,all = T)}
    
}
    return(data_all)}

#单个变量做成list，name_vector是向量
name_list <- function(name_vector=name_vector){
    tem_list <- list()
    for(i in name_vector){
        tem_list[[i]] <- i
    }
    return(tem_list)
}

#批量的多/单因素逻辑回归，data是data(需要先进行numeric)，signle_list是需要批量归一化的彬良列表，
#adj_immo固定校正的变量，比如年龄和性别，group_index是分组变量（为1或者0），
#最后输出的是一个word 的文件，以及tbl的list
#批量logistic函数
adj_log_glm <- function(data=data,
                        signle_list=signle_list,
                        adj_immo=adj_immo,
                        group_index=group_index,
                        outfile=outfile){
    glm_list <- list()
    #如果是空的，就是单因素，如果不空，就是多因素
    if(is.null(adj_immo)){
        for(i in names(signle_list)){
            tem_name_index <- c(i,group_index)
            formul_index <- formula(paste0(group_index,"~."))
            ad_mol <- glm(formul_index,data[,tem_name_index],family=binomial)
            tbl_glm_result <- tbl_regression(ad_mol,exponentiate = T)
            glm_list[[i]] <- tbl_glm_result
        } 
    }else{
    for(i in names(signle_list)){
        tem_name_index <- c(i,adj_immo,group_index)
        formul_index <- formula(paste0(group_index,"~."))
        ad_mol <- glm(formul_index,data[,tem_name_index],family=binomial)
        tbl_glm_result <- tbl_regression(ad_mol,exponentiate = T)
        glm_list[[i]] <- tbl_glm_result
        
        tbl_merge(tbls = glm_list,tab_spanner = names(glm_list))|> 
            as_flex_table() |>  save_as_docx(path = outfile)
    }}
    return(glm_list)
}

#data包含了数据，human_ID是唯一识别号的列名，
#diag_col是包含相关诊断的列名字，diag_disease是我们想诊断的疾病（比如高血压，糖尿病，高脂血症等）
diag_disease_fun <- function(data=data,
                         human_ID=human_ID,
                         diag_col=diag_col,
                         diag_disease=diag_disease){
  #首先生成一个human_ID的矩阵
  data <- as.data.frame(data)
  
  diag_dataframe <- data.frame(unique(data[,human_ID]),diag=NA)
  for(i in 1:nrow(diag_dataframe)){
    patient_id <- diag_dataframe[i,1]
    
    patient_all_diag <- paste(data[data[,human_ID]==patient_id,diag_col],collapse = "_")
    
    diag_dataframe[i,2] <- patient_all_diag
  }
  
  for(j in diag_disease){
    diag_dataframe[,j] <- grepl(paste0(".*",j,".*"),diag_dataframe[,2])
    
  }
  colnames(diag_dataframe)[1] <- human_ID
  diag_dataframe[,human_ID] <- as.character(diag_dataframe[,human_ID])
  return(diag_dataframe)
}

#提取特殊的化验或者药物用量
#fun可以是max,min,sum，但是不可以是mean,num代表整理前多少个药物或者化验
extrade_fun_lab <- function(data=data,
                            human_ID = human_ID,
                            lab_col=lab_col,
                            lab_result=lab_result,fun=fun){
  #首先做出矩阵
  data <- as.data.frame(data)
  
  data <- data[!is.na(data[,lab_col]),]
  data <- data[!is.na(data[,lab_result]),]
  data <- data[!is.na(data[,human_ID]),]
  
  patient_dataframe <- data.frame(unique(data[,human_ID]))
  
  unique_herb <- as.vector(unique(data[,lab_col]))

  for (i in unique_herb) {
    patient_dataframe[,i] <- NA
    
    for(j in 1:nrow(patient_dataframe)){
      
      patien_id <- patient_dataframe[,1][j]
      
      if(!is.na(patien_id)){
        
      tem_herb_data <- data[data[,human_ID]==patien_id & data[,lab_col]==i,c(human_ID,lab_col,lab_result)]
      
      tem_herb_data <-na.omit(tem_herb_data)
      tem_herb_data <- tem_herb_data[,lab_result]
      
      if(length(tem_herb_data)>0){
        
        tem_pati_herb <- fun(na.omit(tem_herb_data))
        patient_dataframe[j,i] <- tem_pati_herb}

      }
    }
  }
  return(patient_dataframe)
}


##寻找出最小的检验,或者药物时间,药物的开出时间小于化验时间
min_time <- function(data_bef=data_bef,
                     data_aft=data_aft,
                     human_ID_bef=human_ID_bef,
                     time_ID_bef=time_ID_bef,
                     human_ID_aft=human_ID_aft,
                     time_ID_aft=time_ID_aft){
  
  #data_bef <- data_bef[!is.na(data_bef[,human_ID_bef])]
  #data_bef <- data_bef[!is.na(data_bef[,time_ID_bef])]
  
  #data_aft <- data_aft[!is.na(data_aft[,human_ID_aft])]
  #data_aft <- data_aft[!is.na(data_aft[,time_ID_aft])]
  
  
  inter_human_ID <- data.frame(intersect(data_bef[,human_ID_bef],data_aft[,human_ID_aft]))
  
  wright_time <- c()
  i <- inter_human_ID[1,1]
  for (i in inter_human_ID[,1]) {
    min_time_bef <- data_bef[data_bef[,human_ID_bef] == i,time_ID_bef] %>% arrange()
    min_time_bef1 <- min_time_bef[nrow(min_time_bef),1]
    
    min_time_aft <- data_aft[data_aft[,human_ID_aft] == i,time_ID_aft] %>% arrange()
    min_time_aft1 <- min_time_aft[nrow(min_time_aft),1]
    
    if (min_time_bef1 < min_time_aft1) {
      wright_time <- cbind(wright_time,i)
    }
    
  }
  return(wright_time)
  
}


#min筛选最小值,
unique_ID <- "就诊标识医渡云计算"
min_ID <- "草药处方开立日期时间"

min_time1 <- function(data=data,
                     unique_ID=unique_ID,
                     min_ID=min_ID){
  
  data <- na.omit(data)

  data <- data[order(as.vector(data[,min_ID])[[1]]),]
  unique_IDD <- as.data.frame(unique(data[,unique_ID]))
  
  i <- unique_IDD[1,1]
  
  time_sum <- c()
  for(i in unique_IDD[,1]){
    
    tem_time_data <- data[data[,unique_ID]==i,][1,]
    time_sum <- rbind(tem_time_data,time_sum)
    
  }
  return(time_sum)
}


