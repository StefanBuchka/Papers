# On the optimistic performance evaluationof newly introduced bioinformatic methods
#The code to generate the quantitative results from the Commentary by Buchka et al.
#This code: Analysis 1: level of the paper
#Latest Date: 31.03.2021

rm(list=ls())

library(dplyr)
library(purrr)

ds <- openxlsx::read.xlsx("data_set_of_extracted_data_Buchka_et_al.xlsx")
#each row represents a method and a rank, and the details of from which substudy the rank came from

dim(ds) #[1] 1131   17

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
rm(paps)
########

#Creating a variable that defines each substudy: a substudy is defined by the paper, the metric of comparison and the datasets used
ds <- ds %>% mutate(substudy = paste(paper,metric,dataset_names,sep="_"))

#for convenience: create a substudy indicator, instead of using complete name
ss    <- unique(ds$substudy)
ss_df <- data.frame(substudy=ss,ss_identifier=paste0("s",1:length(ss)))

#replacing the "substudy" name with the identifier
ds <- inner_join(ds,ss_df) %>% select(-substudy) %>% rename(substudy=ss_identifier)

length(unique(ds$substudy)) #[1] 191

#eliminating column of no use
ds <- select(ds,-metric,-dataset_names,-nb_datasets_substudy,-nb_metrics_substudy,-comment,-method_in_paper_as_named,-name_of_method_in_paper)


###
rm(ss)
###

#qc: checking which columns have missing values
sapply(ds,function(x) sum(is.na(x))) %>% keep(.,~.x>0)
#introduced_method 
#548

########
#Analysis 1: only looking at those substudies that involve all methods of the paper 

#eliminating those comparisons where the total number of methods examined < total number of methods in the paper
ds <- filter(ds,nb_methods_paper==nb_methods_substudy)

#if a method appears more than once in one paper we take its average rank
ds <- ds %>% 
  group_by(paper,introduced_method,method,nb_methods_substudy) %>%
  summarise(mean_rank = mean(rank)) %>% ungroup
dim(ds) #[1] 131   5


###########################

#creates vector with all methods that are introduced in a paper
not_neutral_methods <- filter(ds, !is.na(introduced_method)) %>% .$introduced_method %>% unique
length(not_neutral_methods) #[1] 11

#create a list of the methods to which each of these is compared.   
combination_calc <- lapply(not_neutral_methods,function(x) {
  filter(ds,introduced_method==x, method!=introduced_method) %>% 
      .$method   %>%  unique  })

######################

#a function to extract papers for a given non-neutral method and a given other method
extract_pair <- function(nn_meth,other_meth,dframe=ds) {
  
  #filter the papers to only those where the method is a relevant method
  dframe <- filter(dframe,method %in% c(nn_meth,other_meth))
  
  #creates vector of all papers in which the both desired methods are covered.
  relevant_papers <- dframe %>%
    group_by(paper)         %>%
    distinct(method)        %>%
    filter(.,length(method) == 2)
  
  dframe <- filter(dframe,paper %in% relevant_papers$paper) %>% mutate(comparison=paste(nn_meth,other_meth,sep="_"))
  
  #marking with 0 the methods for which the non-neutral method is that method ranked in the row 
  #ie, it may be a neutral comparison, but the newer method is marked with 0
  dframe$new_or_old_method <- "old"
  dframe$new_or_old_method[which(dframe$method == nn_meth)]  <- "new"
  
  #where papers are other than those introducing the newer method, eliminate "neutral", ie neither method is being introduced
  dframe$introduced_method[which(dframe$introduced_method != nn_meth)] <- NA
  
  return(as.data.frame(dframe,stringsAsFactors=F))
  
}

#for each non-neutral method, get a dataframe of papers for each other method it is compared to 
filtered <- mapply(function(nn_meth,ss) {
  
  #getting the dataframe the non-neutral method and for each of the secondary/older methods
  all_combs <- lapply(ss,extract_pair,nn_meth=nn_meth) %>% do.call(rbind,.)
  return(all_combs)
  
},nn_meth=not_neutral_methods,ss=combination_calc,SIMPLIFY=F) %>% do.call(rbind,.)


#quality control: which columns have NAs
sapply(filtered,function(x) sum(is.na(x))) %>% keep(.,function(x) x>0)
#introduced_method 
#154

#that for every comparison, there is one old and one new method:  here at least that they balance
summary(as.factor(filtered$new_or_old_method))
#new old 
#138 138 

########################
#qc: 
#checking that the introduced method is always new
filter(filtered, introduced_method==method) %>% .$new_or_old_method %>% unique #[1] new
#########################


#selects all comparisons which are compared in two or more papers
filtered <- group_by(filtered,comparison) %>%
  filter(.,length(unique(paper))>1)       %>% 
  ungroup
dim(filtered) #[1] 192   7

#################################
#qc
sapply(unique(filtered$comparison),function(comp) {
  dframe <- filter(filtered,comparison==comp)
  mthds  <- unique(dframe$method)
  qc_1   <- length(mthds)==2                                             #always two methods                
  qc_2   <- sum(dframe$method==mthds[1]) == sum(dframe$method==mthds[2]) #every paper has a rank of both methods
  dframe <- filter(dframe,!is.na(introduced_method))
  qc_3   <- (nrow(dframe) == 2)                                          #both ranks from the paper introducing the newer method
  return(qc_1 & qc_2 & qc_3)
}) %>% all %>% stopifnot
#################################

#Selects all comparisons within a paper in which both desired methods are included: 
#actually unnecessary given the above, but here it is anyway
filtered <- filtered                                 %>%
  group_by(comparison,paper,nb_methods_substudy) %>%
  filter(.,length(method) == 2)                      %>%
  arrange(comparison)                                %>%
  ungroup()

dim(filtered) ##[1] 192   7


##################################
#Goal:
#for each comparison, we get 4 pieces of information:
#1) verify that there is only one paper introducing the new method
#2) whether the newer method is ranked better in the paper introducing it
#3) the number of neutral papers with the comparison
#4) the number of neutral papers with the newer method better than the older

comparisons <- unique(filtered$comparison)
length(comparisons) #[1] 19

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

#given a paper (and one comparison), and whether the paper is introducing a new method,
#we check that the dataframe contains two sensible ranks
paper_qc <- function(dframe,introducing_paper=T) {
  
  qc_1 <- (nrow(dframe)==2)                                                               
  qc_2 <- dframe$introduced_method %>% unique %>% length %>% "=="(1)                      #one introduced method (could be NA)
  qc_3 <- dframe$method            %>% unique %>% length %>% "=="(2)                      #two methods total 
  qc_4 <- setequal(c("new","old"),dframe$new_or_old_method)                               #one method new, one old 
  if (introducing_paper)  { qc_5 <- filter(dframe,introduced_method==method) %>% .$new_or_old_method %>% "=="("new") } #in an introducing paper, that the introduced method is "new"
  if (!introducing_paper) { qc_5 <- filter(dframe,introduced_method==method) %>% nrow(.) %>% "=="(0)                 } #in a non-introducing paper, no method is being introduced
  
  return(qc_1 & qc_2 & qc_3 & qc_4 & qc_5)
}


#################
#From goals above: 

###
#1) verify that there is only one paper introducing the new method AND
#2) if the paper introducing the new method ranks it as better

introducing_papers <- lapply(comparisons,function(comp) {
  
  dframe <- filter(filtered,comparison==comp,!is.na(introduced_method))
  dframe %>% paper_qc %>% stopifnot
  
  #the number of papers
  num_papers <- length(unique(dframe$paper))
  
  #the rankings
  rank_new <- dframe$mean_rank[dframe$new_or_old_method=="new"]
  rank_old <- dframe$mean_rank[dframe$new_or_old_method=="old"]
  
  new_better <- rank_output(rank_new = rank_new,rank_old = rank_old)
  
  return(c(num_papers,new_better))    
})

###
#3) the number of neutral papers with the comparison AND
#4) the number of neutral papers with the newer method better than the older

neutral_papers <- lapply(comparisons,function(comp) {

  #just the neutral papers for the given comparison  
  dframe <- filter(filtered,comparison==comp,is.na(introduced_method))

  #each individual paper
  indiv_papers <- unique(dframe$paper)

  #running qc for each 
  sapply(indiv_papers,function(pap) {
    dfr <- filter(dframe,paper==pap)  
    dfr %>% paper_qc(introducing_paper = F)    
  }) %>% all %>% stopifnot

  num_papers <- length(indiv_papers)

  #determining if the newer method for each paper is superior
  superiority <- sapply(indiv_papers,function(pap) {
    dfr      <- filter(dframe,paper==pap)  
    rank_new <- dfr$mean_rank[dfr$new_or_old_method=="new"]
    rank_old <- dfr$mean_rank[dfr$new_or_old_method=="old"]
    rank_output(rank_new = rank_new,rank_old = rank_old)
  }) 

  return(c(num_papers,sum(superiority)))

})

#final output
comparisons_introducing_papers <- data.frame(comparisons,n_papers=map_dbl(introducing_papers,1),
                                                 n_new_better=map_dbl(introducing_papers,2),stringsAsFactors = F)
comparisons_neutral_papers     <- data.frame(comparisons,n_papers=map_dbl(neutral_papers,1),
                                                 n_new_better=map_dbl(neutral_papers,2),stringsAsFactors = F)

comparisons_papers <- list(introducing = comparisons_introducing_papers,neutral = comparisons_neutral_papers)
comparisons_papers <- lapply(comparisons_papers,function(x) mutate(x,perc_new_better=n_new_better/n_papers))

######
rm(comparisons_introducing_papers,comparisons_neutral_papers)
#####

#Results: 

#total number of comparison pairs
length(comparisons) #[1] 19

#total number of comparisons by type
sapply(comparisons_papers,function(x) sum(x$n_papers))
#introducing     neutral 
#19              77 

#number of times the new method is better
sapply(comparisons_papers,function(x) sum(x$n_new_better))
#introducing     neutral 
#18.0            49.5 

#the percent the newer paper is better
sapply(comparisons_papers,function(x) sum(x$n_new_better)/sum(x$n_papers))
#introducing    neutral 
#0.9473684      0.6428571 

#save the results for the plots
save(comparisons_papers,file="comparisons_papers.RData")






