# On the optimistic performance evaluationof newly introduced bioinformatic methods
#The code to generate the quantitative results from the Commentary by Buchka et al.
#This code: Analysis 2: level of the substudy
#Latest Date: 31.03.2021

rm(list=ls())

library(dplyr)
library(purrr)

ds <- openxlsx::read.xlsx("data_set_of_extracted_data_Buchka_et_al.xlsx")
#each row represents a method and a rank, and the details of from which substudy the rank came from

dim(ds) #[1] 1131   18

#qc: seeing if there are NAs
sapply(ds,function(x) sum(is.na(x))) %>% keep(.,~.x>0)
#introduced_method method_in_paper_as_named  Name_of_method_in_paper        comment 
#548                      360                      771                      506 

#for convenience: create a paper indicator, instead of using complete name
paps    <- unique(ds$paper)
paps_df <- data.frame(paper=paps,pap_identifier=paste0("p",1:length(paps)))

#replacing the "paper" name with the identifier
ds <- inner_join(ds,paps_df) %>% select(-paper) %>% rename(paper=pap_identifier)

length(unique(ds$paper)) #[1] 27

length(unique(ds$introduced_method)) #[1] 16
length(unique(ds$method))            #[1] 97

########
rm(paps,paps_df)
########


#Creating a variable that defines each substudy: a substudy is defined by the paper, the metric of comparison and the datasets used
ds <- mutate(ds,substudy = paste(paper,metric,dataset_names,sep="_"))                                                  


#for convenience: create a substudy indicator, instead of using complete name
ss    <- unique(ds$substudy)
ss_df <- data.frame(substudy=ss,ss_identifier=paste0("s",1:length(ss)))

#replacing the "substudy" name with the identifier
ds <- inner_join(ds,ss_df) %>% select(-substudy) %>% rename(substudy=ss_identifier)

length(unique(ds$substudy)) #[1] 191

#eliminating column of no use
ds <- select(ds,-metric,-dataset_names,-nb_datasets_substudy,-nb_metrics_substudy,-comment,-method_in_paper_as_named,-name_of_method_in_paper)



###
rm(ss,ss_df)
###

#qc: checking which columns have missing values
sapply(ds,function(x) sum(is.na(x))) %>% keep(.,~.x>0)
#introduced_method 
#548


###########################
###########################

#creates vector with all methods that are introduced in a paper
not_neutral_methods <- filter(ds, !is.na(introduced_method)) %>% .$introduced_method %>% unique
length(not_neutral_methods) #[1] 15

#create a list of the methods to which each of these is compared.   
combination_calc <- lapply(not_neutral_methods,function(x) {
  filter(ds,introduced_method==x, method!=introduced_method) %>% 
      .$method   %>%  unique  })

names(combination_calc) <- not_neutral_methods

######################

#a function to extract substudies for a given non-neutral method and a given other method
extract_pair <- function(nn_meth,other_meth,dframe=ds) {
  
  #filter the substudies to only those where the method is a relevant method
  dframe <- filter(dframe,method %in% c(nn_meth,other_meth))
  
  #creates vector of all substudies in which the both desired methods are covered.
  relevant_substudies <- dframe %>%
    group_by(substudy)         %>%
    distinct(method)        %>%
    filter(.,length(method) == 2)
  
  dframe <- filter(dframe,substudy %in% relevant_substudies$substudy) %>% mutate(comparison=paste(nn_meth,other_meth,sep="_"))
  
  if (nrow(dframe)==0) return(NULL)
  
  #marking with 0 the methods for which the non-neutral method is that method ranked in the row 
  #ie, it may be a neutral comparison, but the newer method is marked with 0
  dframe$new_or_old_method <- "old"
  dframe$new_or_old_method[which(dframe$method == nn_meth)]  <- "new"
  
  #where substudies are other than those introducing the newer method, eliminate "neutral", ie neither method is being introduced
  dframe$introduced_method[which(dframe$introduced_method != nn_meth)] <- NA
  
  return(as.data.frame(dframe,stringsAsFactors=F))
  
}

#for each non-neutral method, get a dataframe of substudies for each other method it is compared to 
filtered <- mapply(function(nn_meth,ss) {
  
  #getting the dataframe the non-neutral method and for each of the secondary/older methods
  all_combs <- lapply(ss,extract_pair,nn_meth=nn_meth) %>% do.call(rbind,.)
  return(all_combs)
  
},nn_meth=not_neutral_methods,ss=combination_calc,SIMPLIFY=F) %>% do.call(rbind,.)


#quality control
sapply(filtered,function(x) sum(is.na(x))) %>% keep(.,function(x) x>0)
#introduced_method 
#1342
summary(as.factor(filtered$new_or_old_method)) #that the new methods and old methods balance
#new  old 
#1116 1116

########################
#qc: 
#checking that the introduced method is always new
filter(filtered, introduced_method==method) %>% .$new_or_old_method %>% unique #[1] new
#########################

#selects all comparisons which are compared in two or more substudies
filtered <- group_by(filtered,comparison) %>%
  filter(.,sum(is.na(introduced_method))>0,sum(!is.na(introduced_method))>0) %>% 
  ungroup
dim(filtered) #[1] 1694    14

#################################
#qc
sapply(unique(filtered$comparison),function(comp) {

  dframe <- filter(filtered,comparison==comp)
  mthds  <- unique(dframe$method)
  qc_1   <- length(mthds)==2                                             #always two methods                
  qc_2   <- sum(dframe$method==mthds[1]) == sum(dframe$method==mthds[2]) #every substudy has a rank of both methods
  dframe <- filter(dframe,!is.na(introduced_method))
  qc_3   <- (nrow(dframe) %%  2) == 0                                    #both ranks from the substudy introducing the newer method
  return(qc_1 & qc_2 & qc_3)
}) %>% all %>% stopifnot
#################################

#Selects all comparisons within a substudy in which both desired methods are included: 
#actually unnecessary given the above, but here it is anyway
filtered <- filtered                                 %>%
  group_by(comparison,substudy,nb_methods_substudy)  %>%
  filter(.,length(method) == 2)                      %>%
  arrange(comparison)                                %>%
  ungroup()

dim(filtered) ##[1] 1694   13

##################################

# Exclusion of pair "funnorm_funnormw/noob"
# This is required, because both methods (funnorm and funnormw/noob) were introduced in the same paper. A comparison of pairs
# like this does not make sense. 

filtered <- filtered[-which(filtered$comparison == "funnorm_funnormw/noob"),]
dim(filtered) #[1] 1662   13

##################################
#Goal:
#for each comparison, we get 4 pieces of information:
#1) the number of non-neutral substudies with the comparison
#2) the number of non-neutral substudies with the newer method better than the older
#3) the number of neutral substudies with the comparison
#4) the number of neutral substudies with the newer method better than the older

comparisons <- unique(filtered$comparison)
length(comparisons) #[1] 28

#qc: no comparison is the reverse of any other
rev_comp_a <- strsplit(x=comparisons,split="_") %>% map(1) %>% unlist
rev_comp_b <- strsplit(x=comparisons,split="_") %>% map(2) %>% unlist
rev_comp   <- paste(rev_comp_b,rev_comp_a,sep="_")
intersect(rev_comp,comparisons) %>% length %>% "=="(.,0) %>% stopifnot
rm(rev_comp,rev_comp_a,rev_comp_b)

####################
#important functions

#get an output based on the rank of a new method and an older method
rank_output <- function(rank_new,rank_old) {
  if (rank_new < rank_old)  return(1)
  if (rank_new > rank_old)  return(0)
  if (rank_new == rank_old) return(0.5)
  return(NA)
}

#here: given a substudy (and implicitly one comparison), and whether the substudy is introducing a new method,
#we check that the dataframe contains two sensible ranks
substudy_qc <- function(dframe,introducing_substudy=T) {
  
  qc_1 <- (nrow(dframe)==2)                                                               
  qc_2 <- dframe$introduced_method %>% unique %>% length %>% "=="(1)                  #one introduced method (could be NA)
  qc_3 <- dframe$method            %>% unique %>% length %>% "=="(2)                  #always two methods    
  qc_4 <- setequal(c("new","old"),dframe$new_or_old_method)                           #one method new, one old      
  if (introducing_substudy)  { qc_5 <- filter(dframe,introduced_method==method) %>% .$new_or_old_method %>% "=="("new") } #in an introducing paper, that the introduced method is "new"
  if (!introducing_substudy) { qc_5 <- filter(dframe,introduced_method==method) %>% nrow(.) %>% "=="(0)                 } #in a non-introducing paper, no method is being introduced
  
  return(qc_1 & qc_2 & qc_3 & qc_4 & qc_5)
}



#From goals above: 

###
#1) verify that there is only one substudy introducing the new method AND
#2) if the substudy introducing the new method ranks it as better

introducing_substudies <- lapply(comparisons,function(comp) {
  
  #just the neutral substudies for the given comparison  
  dframe <- filter(filtered,comparison==comp,!is.na(introduced_method))
  
  #each individual substudy
  indiv_substudies <- unique(dframe$substudy)
  
  #running qc for each 
  sapply(indiv_substudies,function(pap) {
    dfr <- filter(dframe,substudy==pap)  
    dfr %>% substudy_qc(introducing_substudy = T)    
  }) %>% all %>% stopifnot
  
  num_substudies <- length(indiv_substudies)
  
  #determining if the newer method for each substudy is superior
  superiority <- sapply(indiv_substudies,function(pap) {
    dfr      <- filter(dframe,substudy==pap)  
    rank_new <- dfr$rank[dfr$new_or_old_method=="new"]
    rank_old <- dfr$rank[dfr$new_or_old_method=="old"]
    rank_output(rank_new = rank_new,rank_old = rank_old)
  }) 
  
  return(c(num_substudies,sum(superiority)))
  
})




###
#3) the number of neutral substudies with the comparison AND
#4) the number of neutral substudies with the newer method better than the older

neutral_substudies <- lapply(comparisons,function(comp) {
    
    #just the neutral substudies for the given comparison  
    dframe <- filter(filtered,comparison==comp,is.na(introduced_method))
    
    #each individual substudy
    indiv_substudies <- unique(dframe$substudy)
    
    #running qc for each 
    sapply(indiv_substudies,function(pap) {
      dfr <- filter(dframe,substudy==pap)  
      dfr %>% substudy_qc(introducing_substudy = F)    
    }) %>% all %>% stopifnot
    
    num_substudies <- length(indiv_substudies)
    
    #determining if the newer method for each substudies is superior
    superiority <- sapply(indiv_substudies,function(pap) {
      dfr      <- filter(dframe,substudy==pap)  
      rank_new <- dfr$rank[dfr$new_or_old_method=="new"]
      rank_old <- dfr$rank[dfr$new_or_old_method=="old"]
      rank_output(rank_new = rank_new,rank_old = rank_old)
    }) 
    
    return(c(num_substudies,sum(superiority)))
    
})


#final output
comparisons_introducing_substudies <- data.frame(comparisons,n_substudies=map_dbl(introducing_substudies,1),
                                                 n_new_better=map_dbl(introducing_substudies,2),stringsAsFactors = F)
comparisons_neutral_substudies     <- data.frame(comparisons,n_substudies=map_dbl(neutral_substudies,1),
                                                 n_new_better=map_dbl(neutral_substudies,2),stringsAsFactors = F)

comparisons_substudies <- list(introducing = comparisons_introducing_substudies,neutral = comparisons_neutral_substudies)
comparisons_substudies <- lapply(comparisons_substudies,function(x) mutate(x,perc_new_better=n_new_better/n_substudies))

######
rm(comparisons_introducing_substudies,comparisons_neutral_substudies)
#####

#Results: 

#total number of comparison pairs
length(comparisons) #[1] 28

#total numbers of substudies in introducing paper and neutral paper
sapply(comparisons_substudies,function(x) sum(x$n_substudies))
# introducing     neutral 
# 164             667 

#substudies where the new method is ranked as better
sapply(comparisons_substudies,function(x) sum(x$n_new_better))
# introducing     neutral 
# 136.5           408.5 

#percentage where the new one is ranked as better
sapply(comparisons_substudies,function(x) sum(x$n_new_better)/sum(x$n_substudies))
# introducing     neutral 
# 0.8323171       0.6124438  


#save the results for the plots
save(comparisons_substudies,file="comparisons_substudies.RData")

