

#_________________Data sets creation of ARR and MRI
#_________________Stefan Buchka
#_________________24.05.2022

library(sas7bdat)
library(tidyverse)
library(haven)

#_______1.) CAMMS223
#####

#_____CREATION OF BASLIBE DATASET

file_base <- "E:/SDT-RP-11223/files_SAN-CAMMS223/Files/Datasets_CAMMS223/adsl.sas7bdat"
camms223_base <- haven::read_sas(data_file = file_base)

base <- data.frame(   study = camms223_base$STUDYID,
                      id = camms223_base$SUBJID,
                      arm = camms223_base$TRTA,
                      age = camms223_base$AGE,
                      sex = camms223_base$SEX,
                      race = camms223_base$RACE,
                      start_date = camms223_base$RFSTDT,
                      end_date = camms223_base$RFENDT,
                      analysis = camms223_base$FAFL, # full analysis population flag
                      completion = NA,
                      region = camms223_base$REGIONDI,
                      dis_duration = camms223_base$TMSIEP * 365.25) %>% 
    mutate(follow_up = as.Date(end_date) - as.Date(start_date)) %>%
    distinct() %>%
  mutate(dis_duration = as.difftime(dis_duration,units = "days")) %>%
  filter(analysis == "Y")


#_____MRI

file_mri <- "E:/SDT-RP-11223/files_SAN-CAMMS223/Files/Datasets_CAMMS223/admr.sas7bdat"
camms223_mri <- haven::read_sas(data_file = file_mri)


mri_temp <- data.frame(study = camms223_mri$STUDYID,
                       id = camms223_mri$SUBJID,
                       arm = camms223_mri$TRTA,
                       para = camms223_mri$PARAM,
                       value = camms223_mri$AVAL,
                       unit = camms223_mri$UNIT,
                       visit = camms223_mri$AVISITN,
                       day = camms223_mri$ANLDY) 


# Take T2 lesion volume
# pivotize data set to wide format

mri_temp <- merge(mri_temp,
                  base,
                  by = c("study","id","arm"),
                  all.y = T) #excludes all patients, which should not be analysed


mri <- mri_temp %>%
  filter(para == "T2 Lesion Load Volume") %>%
  pivot_wider(names_from = para,values_from = value, names_glue = "t2_volume")


#_____ARR  

file_arr <- "E:/SDT-RP-11223/files_SAN-CAMMS223/Files/Datasets_CAMMS223/adrel.sas7bdat"
camms223_arr <- haven::read_sas(data_file = file_arr)

arr_temp <- data.frame(study = camms223_arr$STUDYID,
                       id = camms223_arr$SUBJID,
                       arm = camms223_arr$TRTA,
                       para = camms223_arr$RECAT,
                       date_rel_start = camms223_arr$REDT,
                       date_rel_end = camms223_arr$REENDDTC,
                       relapse = camms223_arr$HADRE,
                       day = camms223_arr$REFLTIME,
                       age = camms223_arr$AGE,
                       sex = camms223_arr$SEX,
                       race = camms223_arr$RACE,
                       start_date = camms223_arr$RFSTDT,
                       end_date = camms223_arr$RFENDT,
                       analysis = camms223_arr$FAFL,
                       region = camms223_arr$REGIONDI) %>%
  distinct()

#____ Number of relapses 2 years befor study start date

pre_stud_rel <- arr_temp %>%
  filter(para == "CLINICAL EPISODE HISTORY") %>%
  filter(analysis == "Y") %>% #exclude non Analysis/ITT population
  mutate(diff = as.Date(start_date) - as.Date(date_rel_start)) %>%
  filter(diff <= 2*365.25) %>%
  select(-diff) %>%
  group_by(id) %>%
  mutate(n_pre_rel = n()) %>%
  select(study,id,arm,age,sex,race,analysis,region,start_date,end_date,n_pre_rel)

# Filter all relapses (relapse == 1) within clinicla trial (para == CLINICAL EPISODE)
# compute follow-up in days
# compute ARR
arr_temp <- arr_temp %>%
  filter(relapse == 1 & para == "CLINICAL EPISODE") %>%
  filter(analysis == "Y") %>% #exclude non Analysis/ITT population
  mutate(follow_up = max(day,na.rm = T)) %>% #Follow-up is up to 5 years, but according to study protocol, relapses are counted until year three only. 
  group_by(id) %>%
  mutate(comulated_rel = n(),
         arr = comulated_rel/(follow_up/365.25)) %>%
  distinct()



arr <- merge(arr_temp,
              pre_stud_rel,
              by = c("study","id","arm","age","sex","race","analysis","region","start_date","end_date"),
             all = T)

arr <- merge(arr,
             base[,-which(names(base) == "follow_up")],
             by = c("study","id","arm","age","sex","race","analysis","region","start_date","end_date"),
             all = T) %>%
  distinct()

#____EDSS

file_edss <- "E:/SDT-RP-11223/files_SAN-CAMMS223/Files/Datasets_CAMMS223/adns.sas7bdat"
camms223_edss <- haven::read_sas(data_file = file_edss)

edss_temp <- data.frame(study = camms223_edss$STUDYID,
                       id = camms223_edss$SUBJID,
                       arm = camms223_edss$TRTA,
                       edss_score = camms223_edss$AVAL,
                       day = camms223_edss$ANLDY,
                       visit = camms223_edss$AVISITN)


edss_temp <- merge(edss_temp,
                  base,
                  by = c("study","id","arm"),
                  all.y = T) #excludes all patients, which should not be analysed


edss <- edss_temp


#_____Disease progression

file_dp <- "E:/SDT-RP-11223/files_SAN-CAMMS223/Files/Datasets_CAMMS223/adms.sas7bdat"
camms223_dp <- haven::read_sas(data_file = file_dp)

dp <- data.frame(study = camms223_dp$STUDYID,
                        id = camms223_dp$SUBJID,
                        arm = camms223_dp$TRTA,
                        para = camms223_dp$PARAM,
                        value = camms223_dp$AVAL,
                        day = camms223_dp$ANLDY,
                        visit = camms223_dp$AVISITN)

dp <- dp %>%
  filter(para %in% c("Average 9 Hole Peg Test","Average 25-Foot Walk")) %>%
  pivot_wider(names_from = para,values_from = value)

names(dp)[c(6,7)] <- c("peg_test","x25_walk")


dp <- merge(dp,
            edss_temp %>% 
              select(study,id,arm,edss_score,day,visit),
            by = c("study","id","arm","day","visit"),
            all.y = T) #excludes all patients, which should not be analysed


rm(camms223_dp)
rm(file_dp)
rm(camms223_mri)
rm(mri_temp)
rm(file_mri)
rm(camms223_arr)
rm(arr_temp)
rm(file_arr)
rm(camms223_edss)
rm(edss_temp)
rm(file_edss)
rm(pre_stud_rel)
rm(base)
rm(file_base)
rm(camms223_base)
#####

#_______2.) CAMMS323
#####

#_____CREATION OF BASLIBE DATASET

file_base <- "E:/SDT-RP-11223/files_SAN-CAMMS323/Files/Datasets_CAMMS323/adsl.sas7bdat"
camms323_base <- haven::read_sas(data_file = file_base)

base <- data.frame(   study = camms323_base$STUDYID,
                      id = camms323_base$SUBJID,
                      arm = camms323_base$TRTA,
                      age = camms323_base$AGE,
                      sex = camms323_base$SEX,
                      race = camms323_base$RACE,
                      start_date = camms323_base$RFSTDT,
                      end_date = camms323_base$RFENDT,
                      analysis = camms323_base$FAFL, # full analysis population flag (small follow-ups are excluded)
                      completion = NA, # No flag found
                      region = camms323_base$REGIONDI,
                      dis_duration = camms323_base$TMSIEP * 365.25) %>% 
  mutate(follow_up = as.Date(end_date) - as.Date(start_date)) %>%
  distinct() %>%
  mutate(dis_duration = as.difftime(dis_duration,units = "days")) %>%
  filter(analysis == "Y")


#_____MRI


file_mri <- "E:/SDT-RP-11223/files_SAN-CAMMS323/Files/Datasets_CAMMS323/admr.sas7bdat"
camms323_mri <- haven::read_sas(data_file = file_mri)

mri_temp <- data.frame(study = camms323_mri$STUDYID,
                       id = camms323_mri$SUBJID,
                       arm = camms323_mri$TRTA,
                       para = camms323_mri$PARAM,
                       value = camms323_mri$AVAL,
                       unit = camms323_mri$UNIT,
                       visit = camms323_mri$AVISITN,
                       day = camms323_mri$ANLDY)


# Take T2 lesion volume, T2 new or enlarged lesion count
# pivotize data set to wide format
mri_temp <- mri_temp %>%
              filter(para %in% c("T2 Lesion Volume","New / Enlarging T2 Lesion Count")) %>%
              pivot_wider(names_from = c(para,unit),values_from = value)


names(mri_temp)[c(6,7)] <- c("t2_new_or_enlarged","t2_volume")
mri_temp <- mri_temp %>%
  mutate(unit = ifelse(!is.na(t2_volume),"cm?",NA))

mri_temp <- merge(mri_temp,
                  base,
                  by = c("study","id","arm"),
                  all.y = T) #excludes all patients, which should not be analysed
  
mri$t2_new_or_enlarged <- NA
mri <- rbind(mri,mri_temp[,names(mri)])

#_____ARR  

file_arr <- "E:/SDT-RP-11223/files_SAN-CAMMS323/Files/Datasets_CAMMS323/adrel.sas7bdat"
camms323_arr <- haven::read_sas(data_file = file_arr)

arr_temp <- data.frame(study = camms323_arr$STUDYID,
                       id = camms323_arr$SUBJID,
                       arm = camms323_arr$TRTA,
                       para = camms323_arr$RECAT,
                       date_rel_start = camms323_arr$REDT,
                       date_rel_end = camms323_arr$REENDDTC,
                       relapse = camms323_arr$HADRE,
                       day = camms323_arr$REFLTIME,
                       age = camms323_arr$AGE,
                       sex = camms323_arr$SEX,
                       race = camms323_arr$RACE,
                       start_date = camms323_arr$RFSTDT,
                       end_date = camms323_arr$RFENDT,
                       analysis = camms323_arr$FAFL, # full analysis population flag
                       region = camms323_arr$REGIONDI)

#Number of relapses 2 years befor study start date
pre_stud_rel <- arr_temp %>%
  filter(para == "CLINICAL EPISODE HISTORY") %>%
  filter(analysis == "Y") %>% #exclude non Analysis/ITT population
  mutate(diff = as.Date(start_date) - as.Date(date_rel_start)) %>%
  filter(diff <= 2*365.25) %>%
  select(-diff) %>%
  group_by(id) %>%
  mutate(n_pre_rel = n()) %>%
  select(study,id,arm,age,sex,race,analysis,region,start_date,end_date,n_pre_rel) %>%
  distinct()


# Filter all relapses (relapse == 1) within clinicla trial (para == CLINICAL EPISODE)
# compute follow-up in days
# compute ARR
arr_temp <- arr_temp %>%
  filter(analysis == "Y") %>% #exclude non Analysis/ITT population
  filter(relapse == 1 & para == "CLINICAL EPISODE") %>%
  mutate(follow_up = as.Date(end_date) - as.Date(start_date)) %>%
  group_by(id) %>%
  mutate(comulated_rel = n(),
         arr = comulated_rel/(as.numeric(follow_up)/365.25)) %>%
  as.data.frame()



arr_temp <- merge(arr_temp,
                  pre_stud_rel,
                  by = c("study","id","arm","age","sex","race","analysis","region","start_date","end_date"),
                  all = T)

arr_temp <- merge(arr_temp[,-which(names(arr_temp) == "follow_up")],
                  base,
                  by = c("study","id","arm","age","sex","race","analysis","region","start_date","end_date"),
                  all = T) %>%
  distinct()


arr <- rbind(arr,arr_temp[,names(arr)])



#____EDSS

file_edss <- "E:/SDT-RP-11223/files_SAN-CAMMS323/Files/Datasets_CAMMS323/adns.sas7bdat"
camms323_edss <- haven::read_sas(data_file = file_edss)

edss_temp <- data.frame(study = camms323_edss$STUDYID,
                        id = camms323_edss$SUBJID,
                        arm = camms323_edss$TRTA,
                        para = camms323_edss$PARAMCD,
                        value = camms323_edss$AVAL,
                        day = camms323_edss$ANLDY,
                        visit = camms323_edss$AVISITN,
                        analysis = camms323_edss$FAFL) %>% # # full analysis population flag 
filter(analysis == "Y")

edss_temp <- edss_temp %>%
  pivot_wider(names_from = para,values_from = value)

names(edss_temp)[7:14] <- c("blad_bow",
                            "brainstm",
                            "cerebell",
                            "cerebral",
                            "edss_score",
                            "pyramidl",
                            "sensory",
                            "visual")

edss <- bind_rows(edss,edss_temp)



#_____Disease progression

file_dp <- "E:/SDT-RP-11223/files_SAN-CAMMS323/Files/Datasets_CAMMS323/adms.sas7bdat"
camms323_dp <- haven::read_sas(data_file = file_dp)

dp_temp <- data.frame(study = camms323_dp$STUDYID,
                 id = camms323_dp$SUBJID,
                 arm = camms323_dp$TRTA,
                 para = camms323_dp$PARAM,
                 value = camms323_dp$AVAL,
                 day = camms323_dp$ANLDY,
                 visit = camms323_dp$AVISITN)

dp_temp <- dp_temp %>%
  filter(para %in% c("Average 9 Hole Peg Test","Average 25-Foot Walk")) %>%
  pivot_wider(names_from = para,values_from = value)

names(dp_temp)[c(6,7)] <- c("x25_walk","peg_test")

dp_temp <- merge(dp_temp,
                 edss_temp %>% 
                   select(study,id,arm,edss_score,pyramidl,day,visit),
                 by = c("study","id","arm","day","visit"),
                 all.y = T)#excludes all patients, which should not be analysed


dp <- bind_rows(dp,dp_temp)

rm(dp_temp)
rm(camms323_dp)
rm(file_dp)
rm(camms323_mri)
rm(mri_temp)
rm(file_mri)
rm(camms323_arr)
rm(arr_temp)
rm(file_arr)
rm(camms323_edss)
rm(edss_temp)
rm(file_edss)
rm(pre_stud_rel)
rm(base)
rm(file_base)
rm(camms323_base)
#####

#_______3.) CAMMS324 
#####

#_____CREATION OF BASLIBE DATASET

file_base <- "E:/SDT-RP-11223/files_SAN-CAMMS32400507/Files/Datasets_CAMMS32400507/adsl.sas7bdat"
camms324_base <- haven::read_sas(data_file = file_base)

base <- data.frame(   study = camms324_base$STUDYID,
                      id = camms324_base$SUBJID,
                      arm = camms324_base$TRTA,
                      age = camms324_base$AGE,
                      sex = camms324_base$SEX,
                      race = camms324_base$RACE,
                      start_date = camms324_base$RFSTDT,
                      end_date = camms324_base$RFENDT,
                      completion = NA,
                      analysis = camms324_base$FAFL, # full analysis population flag (small follow-ups were excluded)
                      region = camms324_base$REGIONDI,
                      dis_duration = camms324_base$TMSIEP * 365.25) %>% 
  mutate(follow_up = as.Date(end_date) - as.Date(start_date)) %>%
  distinct() %>%
  mutate(dis_duration = as.difftime(dis_duration,units = "days")) %>%
  filter(analysis == "Y")

#_____MRI


file_mri <- "E:/SDT-RP-11223/files_SAN-CAMMS32400507/Files/Datasets_CAMMS32400507/admr.sas7bdat"
camms324_mri <- haven::read_sas(data_file = file_mri)

mri_temp <- data.frame(study = camms324_mri$STUDYID,
                       id = camms324_mri$SUBJID,
                       arm = camms324_mri$TRTA,
                       para = camms324_mri$PARAM,
                       value = camms324_mri$AVAL,
                       unit = camms324_mri$UNIT,
                       visit = camms324_mri$AVISITN,
                       day = camms324_mri$ANLDY) %>%
  distinct()


# Take T2 lesion volume, T2 new or enlarged lesion count
# pivotize data set to wide format
mri_temp <- mri_temp %>%
  filter(para %in% c("T2 Lesion Volume","New / Enlarging T2 Lesion Count")) %>%
  pivot_wider(names_from = c(para,unit),values_from = value)


names(mri_temp)[c(6,7)] <- c("t2_new_or_enlarged","t2_volume")
mri_temp <- mri_temp %>%
  mutate(unit = ifelse(!is.na(t2_volume),"cm?",NA))



mri_temp <- merge(mri_temp,
                  base,
                  by = c("study","id","arm"),
                  all.y = T)#excludes all patients, which should not be analysed


mri <- rbind(mri,mri_temp[,names(mri)])

#_____ARR  

file_arr <- "E:/SDT-RP-11223/files_SAN-CAMMS32400507/Files/Datasets_CAMMS32400507/adrel.sas7bdat"
camms324_arr <- haven::read_sas(data_file = file_arr)

arr_temp <- data.frame(study = camms324_arr$STUDYID,
                       id = camms324_arr$SUBJID,
                       arm = camms324_arr$TRTA,
                       para = camms324_arr$RECAT,
                       date_rel_start = camms324_arr$REDT,
                       date_rel_end = camms324_arr$REENDDTC,
                       relapse = camms324_arr$HADRE,
                       day = camms324_arr$REFLTIME,
                       age = camms324_arr$AGE,
                       sex = camms324_arr$SEX,
                       race = camms324_arr$RACE,
                       start_date = camms324_arr$RFSTDT,
                       end_date = camms324_arr$RFENDT,
                       analysis = camms324_arr$FAFL, # full analysis population flag
                       region = camms324_arr$REGIONDI)

#Number of relapses 2 years befor study start date
pre_stud_rel <- arr_temp %>%
  filter(para == "CLINICAL EPISODE HISTORY") %>%
  filter(analysis == "Y") %>% #exclude non Analysis/ITT population
  mutate(diff = as.Date(start_date) - as.Date(date_rel_start)) %>%
  filter(diff <= 2*365.25) %>%
  select(-diff) %>%
  group_by(id) %>%
  mutate(n_pre_rel = n()) %>%
  select(study,id,arm,age,sex,race,analysis,region,start_date,end_date,n_pre_rel) %>%
  distinct()

# Filter all relapses (relapse == 1) within clinicla trial (para == CLINICAL EPISODE)
# compute follow-up in days
# compute ARR
arr_temp <- arr_temp %>%
  filter(relapse == 1 & para == "CLINICAL EPISODE") %>%
  filter(analysis == "Y") %>% #exclude non Analysis/ITT population
  mutate(follow_up = as.Date(end_date) - as.Date(start_date)) %>%
  group_by(id) %>%
  mutate(comulated_rel = n(),
         arr = comulated_rel/(as.numeric(follow_up)/365.25)) %>%
  as.data.frame()



arr_temp <- merge(arr_temp,
                 pre_stud_rel,
                 by = c("study","id","arm","age","sex","race","analysis","region","start_date","end_date"),
                 all = T)

arr_temp <- merge(arr_temp %>% select(-follow_up),
                  base,
                  by = c("study","id","arm","age","sex","race","analysis","region","start_date","end_date"),
                  all = T) %>%
  distinct()

arr <- rbind(arr,arr_temp)

#____EDSS

file_edss <- "E:/SDT-RP-11223/files_SAN-CAMMS32400507/Files/Datasets_CAMMS32400507/adns.sas7bdat"
camms324_edss <- haven::read_sas(data_file = file_edss)

edss_temp <- data.frame(study = camms324_edss$STUDYID,
                        id = camms324_edss$SUBJID,
                        arm = camms324_edss$TRTA,
                        para = camms324_edss$PARAMCD,
                        value = camms324_edss$AVAL,
                        day = camms324_edss$ANLDY,
                        visit = camms324_edss$AVISITN)


edss_temp <- edss_temp %>%
  pivot_wider(names_from = para,values_from = value)

names(edss_temp)[6:13] <- c("blad_bow",
                            "brainstm",
                            "cerebell",
                            "cerebral",
                            "edss_score",
                            "pyramidl",
                            "sensory",
                            "visual")


edss_temp <- merge(edss_temp,
                  base,
                  by = c("study","id","arm"),
                  all.y = T) #excludes all patients, which should not be analysed

edss <- bind_rows(edss,edss_temp)



#_____Disease progression

file_dp <- "E:/SDT-RP-11223/files_SAN-CAMMS32400507/Files/Datasets_CAMMS32400507/adms.sas7bdat"
camms324_dp <- haven::read_sas(data_file = file_dp)

dp_temp <- data.frame(study = camms324_dp$STUDYID,
                      id = camms324_dp$SUBJID,
                      arm = camms324_dp$TRTA,
                      para = camms324_dp$PARAM,
                      value = camms324_dp$AVAL,
                      day = camms324_dp$ANLDY,
                      visit = camms324_dp$AVISITN)

dp_temp <- dp_temp %>%
  filter(para %in% c("Average 9 Hole Peg Test","Average 25-Foot Walk")) %>%
  pivot_wider(names_from = para,values_from = value)

names(dp_temp)[c(6,7)] <- c("peg_test","x25_walk")

dp_temp <- merge(dp_temp,
                 edss_temp %>% 
                   select(study,id,arm,edss_score,pyramidl,day,visit),
                 by = c("study","id","arm","day","visit"),
                 all.y = T) #excludes all patients, which should not be analysed


dp <- bind_rows(dp,dp_temp)



rm(dp_temp)
rm(camms324_dp)
rm(file_dp)
rm(camms324_mri)
rm(mri_temp)
rm(file_mri)
rm(camms324_arr)
rm(arr_temp)
rm(file_arr)
rm(camms324_edss)
rm(edss_temp)
rm(file_edss)
rm(pre_stud_rel)
rm(base)
rm(file_base)
rm(camms324_base)
#####

#_______4.) CFTY720D1201 
#####


#_____BASELINE INFORMATION
# Creation of data set with baseline information

file_base <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D1201/Files/Analysis Ready Datasets/SAS_analysis/a_mshis_mse.sas7bdat"
CFTY720D1201_base <- haven::read_sas(data_file = file_base)

base_temp <- data.frame(study = CFTY720D1201_base$STY1A,
                        id = CFTY720D1201_base$SID1A,
                        arm = CFTY720D1201_base$TGPDSC1A,
                        age = CFTY720D1201_base$age_1n, # age variable looks weird. Only patients aged 28 or 48 included
                        sex = CFTY720D1201_base$SEX1C,
                        race = CFTY720D1201_base$RCE5C,
                        region = CFTY720D1201_base$COU1A,
                        dis_duration = CFTY720D1201_base$DURMS1C*365.25) %>%
  distinct() %>% 
  group_by(id) %>%
  mutate(isdate = ifelse(any(!is.na(dis_duration)),1,0)) %>% 
  filter((!is.na(dis_duration) & isdate == 1)|
           (is.na(dis_duration) & isdate == 0)) %>%
  distinct() %>% 
  arrange(id) %>%
  select(-isdate) %>%
  mutate(dis_duration = as.difftime(dis_duration, units = "days"))


file_base <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D1201/Files/Analysis Ready Datasets/SAS_analysis/a_ident_mse.sas7bdat"
CFTY720D1201_base <- haven::read_sas(data_file = file_base)

base <- data.frame(study = CFTY720D1201_base$STY1A,
                   id = CFTY720D1201_base$SID1A,
                   arm = CFTY720D1201_base$TGPDSC1A,
                   age = CFTY720D1201_base$AGE_1N,
                   sex = CFTY720D1201_base$SEX1C,
                   race = CFTY720D1201_base$RCE5C,
                   start_date = CFTY720D1201_base$FIRSMD1O, # as study start first dose date is defined (frequently used in study protocol)
                   end_date = CFTY720D1201_base$LSTD1O,
                   analysis = CFTY720D1201_base$FAS, #"Full Analysis Set Population": Useless. 
                   completion = CFTY720D1201_base$SBJCMP1C, #Take study completion date instead
                   region = CFTY720D1201_base$COU1A,
                   n_pre_rel = CFTY720D1201_base$NUMRLP2N) %>%
  mutate(analysis = ifelse(analysis == 1,"Y","N"),
         follow_up = as.Date(end_date) - as.Date(start_date) +1) %>%
  filter(completion == 1) %>%
  distinct()

base <- merge(base,
              base_temp,
              by = c("study","id","arm","age","sex","race","region"),
              all.x = T) #excludes all patients, which should not be analysed

# 1 = Caucasian
# 2 = Black
# 3 = Asian
# 4 = Native American
# 5= Pacific Islander
base$race <- "Asian"
base$sex <- ifelse(base$sex == 1,"M","F")


#_____MRI


file_mri <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D1201/Files/Analysis Ready Datasets/SAS_analysis/a_mri_mse.sas7bdat"
CFTY720D1201_mri <- haven::read_sas(data_file = file_mri)

mri_temp <- data.frame(study = CFTY720D1201_mri$STY1A,
                       id = CFTY720D1201_mri$SID1A,
                       arm = CFTY720D1201_mri$TGPDSC1A,
                       t2_new = CFTY720D1201_mri$NUMLES2N,
                       t2_enlarged = CFTY720D1201_mri$NUMLES5N,
                       t2_new_or_enlarged = CFTY720D1201_mri$NUMNNT2,
                       invalid_scan = CFTY720D1201_mri$FSTERPHA,
                       unit = NA,
                       visit = CFTY720D1201_mri$VIS1N,
                       day = CFTY720D1201_mri$DAYNI_1N) %>%
  filter(invalid_scan == 0) %>%
  select(-invalid_scan)


mri_temp$id <- as.factor(mri_temp$id)

mri_temp <- merge(mri_temp,
                  base[,-11], #Exclude number of relapses 2y before study
                  by = c("study","id","arm"),
                  all.y = T) #excludes all patients, which should not be analysed

mri <- dplyr::bind_rows(mri,mri_temp) %>%
  as.data.frame()

#_____ARR  

file_arr <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D1201/Files/Analysis Ready Datasets/SAS_analysis/a_rel_mse.sas7bdat"
CFTY720D1201_arr <- haven::read_sas(data_file = file_arr)

arr_temp <- data.frame(study = CFTY720D1201_arr$STY1A,
                       id = CFTY720D1201_arr$SID1A,
                       arm = CFTY720D1201_arr$TGPDSC1A,
                       date_rel_start = CFTY720D1201_arr$RLPSTT1O,
                       date_rel_end = CFTY720D1201_arr$RLPEND1O,
                       relapse = CFTY720D1201_arr$RLPCFM1C,
                       comulated_rel = CFTY720D1201_arr$RELNUM2, #Number of confirmed relapses
                       day = CFTY720D1201_arr$DAY_1N,
                       follow_up = CFTY720D1201_arr$INSTUDY1, #days in study (non imputed)
                       arr_from_study = CFTY720D1201_arr$ARR1) # ARR of non imputed confirmed relapses


arr_temp <- merge(arr_temp,
                  base,
                  by = c("study","id","arm","follow_up"),
                  all.y = T)#excludes all patients, which should not be analysed

# Filter all relapses (relapse == 1) within clinicla trial (para == CLINICAL EPISODE)
# compute follow-up in days
# compute ARR
arr_temp <- arr_temp %>%
  #filter(analysis == "Y") %>% #exclude non Analysis/ITT population
  mutate(date_rel_start = as.Date(date_rel_start),
         date_rel_end = as.Date(date_rel_end)) %>%
  group_by(id) %>%
  mutate(comulated_rel = sum(relapse, na.rm = T),
         arr = comulated_rel/(follow_up/365.25)) %>%
  filter(relapse > 0 | is.na(relapse))




arr <- arr %>%
  mutate(date_rel_start = as.Date(date_rel_start),
         date_rel_end = as.Date(date_rel_end))


arr <- bind_rows(arr,arr_temp)
arr <- arr %>% select(-para)

#_____EDSS

file_edss <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D1201/Files/Analysis Ready Datasets/SAS_analysis/a_edss_mse.sas7bdat"
CFTY720D1201_edss <- haven::read_sas(data_file = file_edss)

edss_temp <- data.frame(study = CFTY720D1201_edss$STY1A,
                       id = CFTY720D1201_edss$SID1A,
                       arm = CFTY720D1201_edss$TGPDSC1A,
                       edss_score = CFTY720D1201_edss$EXDSC1N,
                       visual = CFTY720D1201_edss$FUNVSU1N,
                       brainstm = CFTY720D1201_edss$FUNBRA1N,
                       pyramidl = CFTY720D1201_edss$FUNPYR1N,
                       cerebell = CFTY720D1201_edss$FUNCER1N,
                       sensory = CFTY720D1201_edss$FUNSEN1N,
                       blad_bow = CFTY720D1201_edss$FUNBOW1N,
                       cerebral = CFTY720D1201_edss$FUNMEN1N,
                       walk = CFTY720D1201_edss$WLKDSN2C,
                       visit = CFTY720D1201_edss$VIS1N,
                       day = CFTY720D1201_edss$DAYNI_1N)


edss_temp$id <- as.factor(edss_temp$id)

edss_temp <- merge(edss_temp,
                  base[,-which(names(base) == "n_pre_rel")],
                  by = c("study","id","arm"),
                  all.y = T) #excludes all patients, which should not be analysed


edss <- dplyr::bind_rows(edss,edss_temp) %>%
  as.data.frame()

#_____Disease Progression

# No information about peg test and 25 ft walk availlable in this study. Another definition for diseas progression will be used. 

dp_temp <- edss %>% select(study,id,arm,day,visit,edss_score,pyramidl)

dp <- bind_rows(dp,dp_temp)

rm(dp_temp)
rm(CFTY720D1201_mri)
rm(mri_temp)
rm(file_mri)
rm(CFTY720D1201_arr)
rm(arr_temp)
rm(file_arr)
rm(CFTY720D1201_edss)
rm(edss_temp)
rm(file_edss)
rm(base)
rm(base_temp)
rm(CFTY720D1201_base)
rm(file_base)
#####

#_______5.) CFTY720D1201E1 
#####

#_____BASELINE INFORMATION
# Creation of data set with baseline information

file_base <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D1201E1/Files/Analysis Ready Datasets/SAS_analysis/a_mshis_mse.sas7bdat"
CFTY720D1201E1_base <- haven::read_sas(data_file = file_base)

base_temp <- data.frame(study = CFTY720D1201E1_base$STY1A,
                        id = CFTY720D1201E1_base$SID1A,
                        arm = CFTY720D1201E1_base$TGPDSC1A,
                        age = CFTY720D1201E1_base$age_1n, # age variable looks weird. Only patients aged 28 or 48 included
                        sex = CFTY720D1201E1_base$SEX1C,
                        race = CFTY720D1201E1_base$RCE5C,
                        region = CFTY720D1201E1_base$COU1A,
                        dis_duration = CFTY720D1201E1_base$DURMS1C*365.25) %>%
  distinct() %>% 
  group_by(id) %>%
  mutate(isdate = ifelse(any(!is.na(dis_duration)),1,0)) %>% 
  filter((!is.na(dis_duration) & isdate == 1)|
           (is.na(dis_duration) & isdate == 0)) %>%
  distinct() %>% 
  arrange(id) %>%
  select(-isdate) %>%
  mutate(dis_duration = as.difftime(dis_duration, units = "days"))


file_base <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D1201E1/Files/Analysis Ready Datasets/SAS_analysis/a_ident_mse.sas7bdat"
CFTY720D1201E1_base <- haven::read_sas(data_file = file_base)

base <- data.frame(study = CFTY720D1201E1_base$STY1A,
                   id = CFTY720D1201E1_base$SID1A,
                   arm = CFTY720D1201E1_base$TGPDSC1A,
                   age = CFTY720D1201E1_base$AGE_1N,
                   sex = CFTY720D1201E1_base$SEX1C,
                   race = CFTY720D1201E1_base$RCE5C,
                   start_date = CFTY720D1201E1_base$FIRSMD1O, # as study start first dose date is defined (frequently used in study protocol)
                   end_date = CFTY720D1201E1_base$LSTSTR1O,
                   analysis = CFTY720D1201E1_base$FAS, #"Full Analysis Set Population": Useless. But Study will be excluded from analysis
                   completion = NA, 
                   region = CFTY720D1201E1_base$COU1A,
                   n_pre_rel = CFTY720D1201E1_base$NUMRLP2N) %>%
  mutate(analysis = ifelse(analysis == 1,"Y","N"),
         follow_up = as.Date(end_date) - as.Date(start_date) +1) %>%
  distinct()

base <- merge(base,
              base_temp,
              by = c("study","id","arm","age","sex","race","region"),
              all = T)


# 1 = Caucasian
# 2 = Black
# 3 = Asian
# 4 = Native American
# 5= Pacific Islander
base$race <- "Asian"
base$sex <- ifelse(base$sex == 1,"M","F")
#base$study <- "CFTY720D1201E"


#_____MRI

file_mri <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D1201E1/Files/Analysis Ready Datasets/SAS_analysis/a_mri_mse.sas7bdat"
CFTY720D1201E1_mri <- haven::read_sas(data_file = file_mri)

mri_temp <- data.frame(study = CFTY720D1201E1_mri$STY1A,
                       id = CFTY720D1201E1_mri$SID1A,
                       arm = CFTY720D1201E1_mri$TGPDSC1A,
                       t2_new = CFTY720D1201E1_mri$NUMLES2N,
                       t2_enlarged = CFTY720D1201E1_mri$NUMLES5N,
                       t2_new_or_enlarged = CFTY720D1201E1_mri$NUMNNT2,
                       invalid_scan = CFTY720D1201E1_mri$FSTERPHA,
                       unit = NA,
                       visit = CFTY720D1201E1_mri$VIS1N,
                       day = CFTY720D1201E1_mri$DAYNI_1N) %>%
  filter(invalid_scan == 0) %>%
  select(-invalid_scan)


mri_temp$id <- as.factor(mri_temp$id)

mri_temp <- merge(mri_temp,
                  base[,-which(names(base) == "n_pre_rel")],
                  by = c("study","id","arm"),
                  all = T)

mri_temp$study <- "CFTY720D1201E"



mri <- dplyr::bind_rows(mri,mri_temp) %>%
  as.data.frame()

#_____ARR  

file_arr <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D1201E1/Files/Analysis Ready Datasets/SAS_analysis/a_rel_mse.sas7bdat"
CFTY720D1201E1_arr <- haven::read_sas(data_file = file_arr)

arr_temp <- data.frame(study = CFTY720D1201E1_arr$STY1A,
                       id = CFTY720D1201E1_arr$SID1A,
                       arm = CFTY720D1201E1_arr$TGPDSC1A,
                       date_rel_start = CFTY720D1201E1_arr$RLPSTT1O,
                       date_rel_end = CFTY720D1201E1_arr$RLPEND1O,
                       relapse = CFTY720D1201E1_arr$RLPCFM1C,
                       comulated_rel = CFTY720D1201E1_arr$RELNUM2, #Number of confirmed relapses
                       day = CFTY720D1201E1_arr$DAY_1N,
                       follow_up = CFTY720D1201E1_arr$INSTUDY1, #days in study (non imputed)
                       arr_from_study = CFTY720D1201E1_arr$ARR1) # ARR of non imputed confirmed relapses


arr_temp <- merge(arr_temp,
                  base,
                  by = c("study","id","arm","follow_up"),
                  all = T)

# Filter all relapses (relapse == 1) within clinicla trial (para == CLINICAL EPISODE)
# compute follow-up in days
# compute ARR
arr_temp <- arr_temp %>%
  filter(analysis == "Y") %>% #exclude non Analysis/ITT population
  mutate(date_rel_start = as.Date(date_rel_start),
         date_rel_end = as.Date(date_rel_end)) %>%
  group_by(id) %>%
  mutate(comulated_rel = sum(relapse, na.rm = T),
         arr = comulated_rel/(follow_up/365.25)) %>%
  distinct() %>%
  filter(relapse > 0 | is.na(relapse))

arr_temp$study <- "CFTY720D1201E"



arr <- bind_rows(arr,arr_temp)

#_____EDSS

file_edss <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D1201E1/Files/Analysis Ready Datasets/SAS_analysis/a_edss_mse.sas7bdat"
CFTY720D1201E1_edss <- haven::read_sas(data_file = file_edss)

edss_temp <- data.frame(study = CFTY720D1201E1_edss$STY1A,
                        id = CFTY720D1201E1_edss$SID1A,
                        arm = CFTY720D1201E1_edss$TGPDSC1A,
                        edss_score = CFTY720D1201E1_edss$EXDSC1N,
                        visual = CFTY720D1201E1_edss$FUNVSU1N,
                        brainstm = CFTY720D1201E1_edss$FUNBRA1N,
                        pyramidl = CFTY720D1201E1_edss$FUNPYR1N,
                        cerebell = CFTY720D1201E1_edss$FUNCER1N,
                        sensory = CFTY720D1201E1_edss$FUNSEN1N,
                        blad_bow = CFTY720D1201E1_edss$FUNBOW1N,
                        cerebral = CFTY720D1201E1_edss$FUNMEN1N,
                        walk = CFTY720D1201E1_edss$WLKDSN2C,
                        visit = CFTY720D1201E1_edss$VIS1N,
                        day = CFTY720D1201E1_edss$DAYNI_1N)


edss_temp$id <- as.factor(edss_temp$id)



edss_temp <- merge(edss_temp,
                  base[,-which(names(base) == "n_pre_rel")],
                  by = c("study","id","arm"),
                  all = T)

edss_temp$study <- "CFTY720D1201E"


edss <- dplyr::bind_rows(edss,edss_temp) %>%
  as.data.frame()


#_____Diseas Progression

# No information about peg tet and 25 walk test availlable in this study. So another definition for diseas progression will be used


dp_temp <- edss_temp %>%
  select(study,id,arm,visit,day,edss_score,pyramidl)


dp <- bind_rows(dp,dp_temp)

rm(dp_temp)
rm(CFTY720D1201E1_mri)
rm(mri_temp)
rm(file_mri)
rm(CFTY720D1201E1_arr)
rm(arr_temp)
rm(file_arr)
rm(CFTY720D1201E1_edss)
rm(edss_temp)
rm(file_edss)
rm(base)
rm(file_base)
rm(CFTY720D1201E1_base)
rm(base_temp)
#####

#_______6.) CFTY720D2301 
#####

#_____BASELINE INFORMATION
# Creation of data set with baseline information

file_base <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2301/Files/Analysis Ready Datasets/SAS_analysis/a_mshis.sas7bdat"
CFTY720D2301_base <- haven::read_sas(data_file = file_base)

base_temp <- data.frame(study = CFTY720D2301_base$STY1A,
                        id = CFTY720D2301_base$STYSID1A,
                        arm = CFTY720D2301_base$TGPDSC1A,
                        age = NA,
                        age_cat = CFTY720D2301_base$age_1n_cat,
                        sex = CFTY720D2301_base$SEX1C,
                        race = CFTY720D2301_base$RCE5C,
                        region = CFTY720D2301_base$COU1A,
                        dis_date = CFTY720D2301_base$FIRSY_1O,
                        n_pre_rel = CFTY720D2301_base$NUMRLP2N) %>%
  distinct() %>% 
  group_by(id) %>%
  mutate(isdate = ifelse(any(!is.na(dis_date)),1,0)) %>% 
  filter((!is.na(dis_date) & isdate == 1)|
           (is.na(dis_date) & isdate == 0)) %>%
  distinct() %>% 
  mutate(isrel = ifelse(any(!is.na(n_pre_rel)),1,0)) %>% 
  filter((!is.na(n_pre_rel) & isrel == 1)|
           (is.na(n_pre_rel) & isrel == 0)) %>%
  arrange(id) %>%
  select(-isdate, -isrel)


file_base <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2301/Files/Analysis Ready Datasets/SAS_analysis/a_ident.sas7bdat"
CFTY720D2301_base <- haven::read_sas(data_file = file_base)

base <- data.frame(study = CFTY720D2301_base$STY1A,
                   id = CFTY720D2301_base$STYSID1A,
                   arm = CFTY720D2301_base$TGPDSC1A,
                   age = NA,
                   age_cat = CFTY720D2301_base$age_1n_cat,
                   sex = CFTY720D2301_base$SEX1C,
                   race = CFTY720D2301_base$RCE5C,
                   start_date = CFTY720D2301_base$FIRSMD1O, # as study start first dose date is defined (frequently used in study protocol)
                   end_date = CFTY720D2301_base$CUTOFF1O,
                   analysis = NA, #not availlable
                   completion = CFTY720D2301_base$SBJCMP1C, #Take study completion date instead
                   follow_up = CFTY720D2301_base$INSTUDY1,
                   region = CFTY720D2301_base$COU1A) %>%
  mutate(analysis = ifelse(analysis == 1,"Y","N"),
         follow_up = base::as.difftime(follow_up, units = "days")) %>%
  filter(completion == 1) %>%
  distinct()



base <- merge(base,
              base_temp,
              by = c("study","id","arm","age","age_cat","sex","race","region"),
              all.x = T) #excludes all patients, which should not be analysed

base <- base %>%
  mutate(dis_duration = as.Date(start_date) - as.Date(dis_date)) %>%
  select(-dis_date)

# 1 = Caucasian
# 2 = Black
# 3 = Asian
# 7 = Native American
# 8 = Pacific islander
# 77 = Unknown
# 88 = Other

base$race <- as.factor(base$race)
levels(base$race) <- c("Caucasian","Black","Asian","Other")
base$race <- as.character(base$race)
base$sex <- ifelse(base$sex == 1,"M","F")



#_____MRI

file_mri <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2301/Files/Analysis Ready Datasets/SAS_analysis/a_mri.sas7bdat"
CFTY720D2301_mri <- haven::read_sas(data_file = file_mri)

mri_temp <- data.frame(study = CFTY720D2301_mri$STY1A,
                       id = CFTY720D2301_mri$STYSID1A,
                       arm = CFTY720D2301_mri$TGPDSC1A,
                       para = CFTY720D2301_mri$MRIRSL1A,
                       unit = CFTY720D2301_mri$RSLUNT1A,
                       value = CFTY720D2301_mri$RSLVAL1N,
                       invalid_scan = CFTY720D2301_mri$FSTERPHA,
                       visit = CFTY720D2301_mri$VIS1N,
                       day = CFTY720D2301_mri$DAYNI_1N) %>%
  filter(invalid_scan == 0) %>%
  select(-invalid_scan)

# Take T2 lesion volume, T2 new or enlarged lesion count
# pivotize data set to wide format
mri_temp <- mri_temp %>%
  filter(para %in% c("T2_Vol","T2New_R")) %>%
  pivot_wider(names_from = c(para,unit),values_from = value) %>%
  select(-T2_Vol_)


names(mri_temp)[c(6,7)] <- c("t2_volume","t2_new_or_enlarged")
mri_temp <- mri_temp %>%
  mutate(unit = ifelse(!is.na(t2_volume),"mm?",NA))


mri_temp$id <- as.factor(mri_temp$id)

mri_temp <- merge(mri_temp,
                  base[,-which(names(base) == "n_pre_rel")],
                  by = c("study","id","arm"),
                  all.y = T) #excludes all patients, which should not be analysed


mri <- dplyr::bind_rows(mri,mri_temp) %>%
  as.data.frame()

#_____ARR  

file_arr <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2301/Files/Analysis Ready Datasets/SAS_analysis/a_rel.sas7bdat"
CFTY720D2301_arr <- haven::read_sas(data_file = file_arr)

arr_temp <- data.frame(study = CFTY720D2301_arr$STY1A,
                       id = CFTY720D2301_arr$STYSID1A,
                       arm = CFTY720D2301_arr$TGPDSC1A,
                       date_rel_start = CFTY720D2301_arr$RLPSTT1O,
                       date_rel_end = CFTY720D2301_arr$RLPEND1O,
                       relapse = CFTY720D2301_arr$RLPCFM1C,
                       comulated_rel = CFTY720D2301_arr$RELNUM2, #Number of confirmed relapses
                       day = CFTY720D2301_arr$DAY_1N,
                       follow_up = (CFTY720D2301_arr$INSTUDY1) , #days in study (non imputed)
                       arr_from_study = CFTY720D2301_arr$ARR1) # ARR of non imputed confirmed relapses

arr_temp <- merge(arr_temp,
                  base[,c("study","id","arm","follow_up")],
                  by = c("study","id","arm"),
                  all.y = T) #excludes all patients, which should not be analysed


arr_temp[which(arr_temp$follow_up.x != arr_temp$follow_up.y),"follow_up.x"] <- arr_temp[which(arr_temp$follow_up.x != arr_temp$follow_up.y),"follow_up.y"] 
names(arr_temp)[which(names(arr_temp) == "follow_up.x")] <- "follow_up"
arr_temp$follow_up.y <- NULL

arr_temp <- unique(arr_temp)

arr_temp <- merge(arr_temp,
                  base,
                  by = c("study","id","arm","follow_up"),
                  all.y = T) #excludes all patients, which should not be analysed

# Filter all relapses (relapse == 1) within clinicla trial (para == CLINICAL EPISODE)
# compute follow-up in days
# compute ARR
arr_temp <- arr_temp %>%
  #filter(analysis == "Y") %>% #exclude non Analysis/ITT population
  mutate(date_rel_start = as.Date(date_rel_start),
         date_rel_end = as.Date(date_rel_end)) %>%
  group_by(id) %>%
  mutate(comulated_rel = sum(relapse, na.rm = T),
         arr = comulated_rel/(follow_up/365.25)) %>%
  distinct() %>%
  filter(relapse > 0 | is.na(relapse))



arr <- bind_rows(arr,arr_temp)


#_____EDSS

file_edss <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2301/Files/Analysis Ready Datasets/SAS_analysis/a_edss.sas7bdat"
CFTY720D2301_edss <- haven::read_sas(data_file = file_edss)

edss_temp <- data.frame(study = CFTY720D2301_edss$STY1A,
                        id = CFTY720D2301_edss$STYSID1A,
                        arm = CFTY720D2301_edss$TGPDSC1A,
                        edss_score = CFTY720D2301_edss$EXDSC1N,
                        visual = CFTY720D2301_edss$FUNVSU1N,
                        brainstm = CFTY720D2301_edss$FUNBRA1N,
                        pyramidl = CFTY720D2301_edss$FUNPYR1N,
                        cerebell = CFTY720D2301_edss$FUNCER1N,
                        sensory = CFTY720D2301_edss$FUNSEN1N,
                        blad_bow = CFTY720D2301_edss$FUNBOW1N,
                        cerebral = CFTY720D2301_edss$FUNMEN1N,
                        walk = CFTY720D2301_edss$WLKDSN2C,
                        visit = CFTY720D2301_edss$VIS1N,
                        day = CFTY720D2301_edss$DAYNI_1N)


edss_temp$id <- as.factor(edss_temp$id)

edss_temp <- merge(edss_temp,
                  base[,-which(names(base) == "n_pre_rel")],
                  by = c("study","id","arm"),
                  all.y = T) #excludes all patients, which should not be analysed

edss <- dplyr::bind_rows(edss,edss_temp) %>%
  as.data.frame()



#_____Disease progression

file_dp <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2301/Files/Analysis Ready Datasets/SAS_analysis/a_msfc.sas7bdat"
CFTY720D2301_dp <- haven::read_sas(data_file = file_dp)

dp_temp <- data.frame(study = CFTY720D2301_dp$STY1A,
                      id = CFTY720D2301_dp$STYSID1A,
                      arm = CFTY720D2301_dp$TGPDSC1A,
                      peg_test = CFTY720D2301_dp$HPT,
                      x25_walk = CFTY720D2301_dp$T25FW,
                      visit = CFTY720D2301_dp$VIS1N,
                      day = CFTY720D2301_dp$DAYNI_1N) %>%
  unique()

dp_temp <- merge(dp_temp,
                 edss_temp %>% 
                   select(study,id,arm,edss_score,pyramidl,day,visit),
                 by = c("study","id","arm","day","visit"),
                 all.y = T) #excludes all patients, which should not be analysed


dp <- bind_rows(dp,dp_temp)



rm(dp_temp)
rm(CFTY720D2301_dp)
rm(file_dp)
rm(CFTY720D2301_mri)
rm(mri_temp)
rm(file_mri)
rm(CFTY720D2301_arr)
rm(arr_temp)
rm(file_arr)
rm(CFTY720D2301_edss)
rm(edss_temp)
rm(file_edss)
rm(file_base)
rm(base)
rm(CFTY720D2301_base)
rm(base_temp)
#####

#_______7.) CFTY720D2302 
#####

#_____BASELINE INFORMATION
# Creation of data set with baseline information

file_base <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2302/Files/Analysis Ready Datasets/SAS_analysis/a_mshis.sas7bdat"
CFTY720D2302_base <- haven::read_sas(data_file = file_base)

base_temp <- data.frame(study = CFTY720D2302_base$STY1A,
                        id = CFTY720D2302_base$STYSID1A,
                        arm = CFTY720D2302_base$TGPDSC1A,
                        age = NA,
                        age_cat = CFTY720D2302_base$age_1n_cat,
                        sex = CFTY720D2302_base$SEX1C,
                        race = CFTY720D2302_base$RCE5C,
                        region = CFTY720D2302_base$COU1A,
                        dis_date = CFTY720D2302_base$FIRSY_1O,
                        n_pre_rel = CFTY720D2302_base$NUMRLP2N) %>%
  distinct() %>% 
  group_by(id) %>%
  mutate(isdate = ifelse(any(!is.na(dis_date)),1,0)) %>% 
  filter((!is.na(dis_date) & isdate == 1)|
           (is.na(dis_date) & isdate == 0)) %>%
  distinct() %>% 
  mutate(isrel = ifelse(any(!is.na(n_pre_rel)),1,0)) %>% 
  filter((!is.na(n_pre_rel) & isrel == 1)|
           (is.na(n_pre_rel) & isrel == 0)) %>%
  arrange(id) %>%
  select(-isdate, -isrel)


file_base <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2302/Files/Analysis Ready Datasets/SAS_analysis/a_ident.sas7bdat"
CFTY720D2302_base <- haven::read_sas(data_file = file_base)

base <- data.frame(study = CFTY720D2302_base$STY1A,
                   id = CFTY720D2302_base$STYSID1A,
                   arm = CFTY720D2302_base$TGPDSC1A,
                   age = NA,
                   age_cat = CFTY720D2302_base$age_1n_cat,
                   sex = CFTY720D2302_base$SEX1C,
                   race = CFTY720D2302_base$RCE5C,
                   start_date = CFTY720D2302_base$FIRSMD1O, # as study start first dose date is defined (frequently used in study protocol)
                   end_date = NA,
                   end_date_calc = CFTY720D2302_base$FIRSMD2O,
                   follow_up_calc1 = CFTY720D2302_base$SBJCMP1C,
                   follow_up_calc2 = CFTY720D2302_base$SBJDCN2C,
                   follow_up_calc3 = CFTY720D2302_base$FIRSMD2D,
                   analysis = NA, #not availlable
                   completion = CFTY720D2302_base$SBJCMP1C, #Take study completion date instead
                   region = CFTY720D2302_base$COU1A) %>%
  filter(completion == 1)

file_base <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2302/Files/Analysis Ready Datasets/SAS_analysis/a_vis.sas7bdat"
CFTY720D2302_base <- haven::read_sas(data_file = file_base)

base2 <- data.frame(study = CFTY720D2302_base$STY1A,
                   id = CFTY720D2302_base$STYSID1A,
                   arm = CFTY720D2302_base$TGPDSC1A,
                   vis10 = CFTY720D2302_base$VIS1O,
                   vis1n = CFTY720D2302_base$VIS1N) %>%
  filter(vis1n == 778) %>%
  select(-vis1n)

base <- merge(base,
              base2,
              by = c("study","id","arm"),
              all.x = T)#excludes all patients, which should not be analysed

rm(base2)

library(lubridate)

base <- base %>%
  group_by(study,id) %>%
  mutate(analysis = ifelse(analysis == 1,"Y","N"),
         end_date = ifelse(follow_up_calc1 == 1 & follow_up_calc2 != 1 & follow_up_calc3 != "",end_date_calc,end_date), # Procedure according to study protocoll. 
         end_date = ifelse(follow_up_calc3 == "",
                       min(c((start_date + lubridate::days(374)),(vis10 + lubridate::days(1))), na.rm = T),
                       end_date),
         end_date = as.Date(as.POSIXct(end_date, origin = "1970-01-01")),
         follow_up = as.Date(end_date) - as.Date(start_date)) %>%
  distinct() %>% 
  select(-contains("calc"), -vis10) %>% as.data.frame()


base <- merge(base,
              base_temp,
              by = c("study","id","arm","age","age_cat","sex","race","region"),
              all.x = T) #excludes all patients, which should not be analysed

base <- base %>%
  mutate(dis_duration = as.Date(start_date) - as.Date(dis_date)) %>%
  select(-dis_date)

# 1 = Caucasian
# 2 = Black
# 3 = Asian
# 7 = Native American
# 8 = Pacific islander
# 77 = Unknown
# 88 = Other

base$race <- as.factor(base$race)
levels(base$race) <- c("Caucasian","Black","Asian","Native American","Other")
base$race <- as.character(base$race)
base$sex <- ifelse(base$sex == 1,"M","F")

#_____MRI

file_mri <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2302/Files/Analysis Ready Datasets/SAS_analysis/a_mri.sas7bdat"
CFTY720D2302_mri <- haven::read_sas(data_file = file_mri)

mri_temp <- data.frame(study = CFTY720D2302_mri$STY1A,
                       id = CFTY720D2302_mri$STYSID1A,
                       arm = CFTY720D2302_mri$TGPDSC1A,
                       para = CFTY720D2302_mri$MRIRSL1A,
                       unit = CFTY720D2302_mri$RSLUNT1A,
                       value = CFTY720D2302_mri$RSLVAL1N,
                       invalid_scan = CFTY720D2302_mri$FSTERPHA,
                       visit = CFTY720D2302_mri$VIS1N,
                       day = CFTY720D2302_mri$DAYNI_1N) %>%
  filter(invalid_scan == 0) %>%
  select(-invalid_scan)


# Take T2 lesion volume, T2 new or enlarged lesion count
# pivotize data set to wide format
mri_temp <- mri_temp %>%
  filter(para %in% c("PD_VOL","T2_NWEN")) %>%
  pivot_wider(names_from = c(para,unit),values_from = value)


names(mri_temp)[c(6,7)] <- c("t2_volume","t2_new_or_enlarged")
mri_temp <- mri_temp %>%
  mutate(unit = ifelse(!is.na(t2_volume),"mm?",NA))


mri_temp$id <- as.factor(mri_temp$id)

mri_temp <- merge(mri_temp,
                  base[,-which(names(base) == "n_pre_rel")],
                  by = c("study","id","arm"),
                  all.y = T)#excludes all patients, which should not be analysed


mri <- dplyr::bind_rows(mri,mri_temp) %>%
  as.data.frame()

#_____ARR  

file_arr <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2302/Files/Analysis Ready Datasets/SAS_analysis/a_rel.sas7bdat"
CFTY720D2302_arr <- haven::read_sas(data_file = file_arr)

arr_temp <- data.frame(study = CFTY720D2302_arr$STY1A,
                       id = CFTY720D2302_arr$STYSID1A,
                       arm = CFTY720D2302_arr$TGPDSC1A,
                       date_rel_start = CFTY720D2302_arr$RLPSTT1O,
                       date_rel_end = CFTY720D2302_arr$RLPEND1O,
                       relapse = CFTY720D2302_arr$RLPCFM1C,
                       comulated_rel = CFTY720D2302_arr$RELNUM2, #Number of confirmed relapses
                       day = CFTY720D2302_arr$DAYNI_1N,
                       follow_up = CFTY720D2302_arr$INSTUDY1, #days in study (non imputed)
                       arr_from_study = CFTY720D2302_arr$ARR1) # ARR of non imputed confirmed relapses


arr_temp <- merge(arr_temp,
                  base,
                  by = c("study","id","arm","follow_up"),
                  all.y = T)#excludes all patients, which should not be analysed


# Filter all relapses (relapse == 1) within clinicla trial (para == CLINICAL EPISODE)
# compute follow-up in days
# compute ARR
arr_temp <- arr_temp %>%
  #filter(analysis == "Y") %>% #exclude non Analysis/ITT population
  mutate(date_rel_start = as.Date(date_rel_start),
         date_rel_end = as.Date(date_rel_end)) %>%
  group_by(id) %>%
  mutate(comulated_rel = sum(relapse, na.rm = T),
         arr = comulated_rel/(follow_up/365.25)) %>%
  distinct() %>%
  filter(relapse > 0 | is.na(relapse)) %>%
  as.data.frame()


arr <- bind_rows(arr,arr_temp)

#_____EDSS

file_edss <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2302/Files/Analysis Ready Datasets/SAS_analysis/a_edss.sas7bdat"
CFTY720D2302_edss <- haven::read_sas(data_file = file_edss)

edss_temp <- data.frame(study = CFTY720D2302_edss$STY1A,
                        id = CFTY720D2302_edss$STYSID1A,
                        arm = CFTY720D2302_edss$TGPDSC1A,
                        edss_score = CFTY720D2302_edss$EXDSC1N,
                        visual = CFTY720D2302_edss$FUNVSU1N,
                        brainstm = CFTY720D2302_edss$FUNBRA1N,
                        pyramidl = CFTY720D2302_edss$FUNPYR1N,
                        cerebell = CFTY720D2302_edss$FUNCER1N,
                        sensory = CFTY720D2302_edss$FUNSEN1N,
                        blad_bow = CFTY720D2302_edss$FUNBOW1N,
                        cerebral = CFTY720D2302_edss$FUNMEN1N,
                        walk = CFTY720D2302_edss$WLKDSN1C,
                        visit = CFTY720D2302_edss$VIS1N,
                        day = CFTY720D2302_edss$DAYNI_1N)


edss_temp$id <- as.factor(edss_temp$id)

edss_temp <- merge(edss_temp,
                  base[,-which(names(base) == "n_pre_rel")],
                  by = c("study","id","arm"),
                  all.y = T)#excludes all patients, which should not be analysed


edss <- dplyr::bind_rows(edss,edss_temp) %>%
  as.data.frame()



#_____Disease progression

file_dp <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2302/Files/Analysis Ready Datasets/SAS_analysis/a_msfc.sas7bdat"
CFTY720D2302_dp <- haven::read_sas(data_file = file_dp)

dp_temp <- data.frame(study = CFTY720D2302_dp$STY1A,
                      id = CFTY720D2302_dp$STYSID1A,
                      arm = CFTY720D2302_dp$TGPDSC1A,
                      peg_test = CFTY720D2302_dp$HPT,
                      x25_walk = CFTY720D2302_dp$T25FW,
                      visit = CFTY720D2302_dp$VIS1N,
                      day = CFTY720D2302_dp$DAYNI_1N) %>%
  unique()


dp_temp <- merge(dp_temp,
                 edss_temp %>% 
                   select(study,id,arm,edss_score,pyramidl,day,visit),
                 by = c("study","id","arm","day","visit"),
                 all.y = T)#excludes all patients, which should not be analysed

dp <- bind_rows(dp,dp_temp)



rm(dp_temp)
rm(CFTY720D2302_dp)
rm(file_dp)
rm(CFTY720D2302_mri)
rm(mri_temp)
rm(file_mri)
rm(CFTY720D2302_arr)
rm(arr_temp)
rm(file_arr)
rm(CFTY720D2302_edss)
rm(edss_temp)
rm(file_edss)
rm(file_base)
rm(base)
rm(CFTY720D2302_base)
rm(base_temp)
#####

#_______8.) CFTY720D2309
#####

#_____BASELINE INFORMATION
# Creation of data set with baseline information

file_base <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2309/Files/Analysis Ready Datasets/SAS_analysis/a_mshis.sas7bdat"
CFTY720D2309_base <- haven::read_sas(data_file = file_base)

base_temp <- data.frame(study = CFTY720D2309_base$STY1A,
                        id = CFTY720D2309_base$STYSID1A,
                        age = NA,
                        age_cat = CFTY720D2309_base$age_1n_cat,
                        sex = CFTY720D2309_base$SEX1C,
                        race = CFTY720D2309_base$RCE5C,
                        region = CFTY720D2309_base$COU1A,
                        dis_date = CFTY720D2309_base$FIRSY_1O,
                        n_pre_rel = CFTY720D2309_base$NUMRLP2N) %>%
  distinct() %>% 
  group_by(id) %>%
  mutate(isdate = ifelse(any(!is.na(dis_date)),1,0)) %>% 
  filter((!is.na(dis_date) & isdate == 1)|
           (is.na(dis_date) & isdate == 0)) %>%
  distinct() %>% 
  mutate(isrel = ifelse(any(!is.na(n_pre_rel)),1,0)) %>% 
  filter((!is.na(n_pre_rel) & isrel == 1)|
           (is.na(n_pre_rel) & isrel == 0)) %>%
  arrange(id) %>%
  select(-isdate, -isrel)


file_base <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2309/Files/Analysis Ready Datasets/SAS_analysis/a_ident.sas7bdat"
CFTY720D2309_base <- haven::read_sas(data_file = file_base)

base <- data.frame(study = CFTY720D2309_base$sty1a,
                   id = CFTY720D2309_base$stysid1a,
                   arm = CFTY720D2309_base$tgpdsc1a,
                   age = NA,
                   age_cat = CFTY720D2309_base$age_1n_cat,
                   sex = CFTY720D2309_base$sex1c,
                   race = CFTY720D2309_base$rce5c,
                   start_date = CFTY720D2309_base$firsmd1o, # as study start first dose date is defined (frequently used in study protocol)
                   end_date = CFTY720D2309_base$cutoff1o,
                   analysis = NA, 
                   completion = CFTY720D2309_base$sbjcmp1c, #Take study completion date instead
                   follow_up = CFTY720D2309_base$instudy1,
                   region = CFTY720D2309_base$cou1a) %>%
  mutate(analysis = ifelse(analysis == 1,"Y","N"),
         follow_up = base::as.difftime(follow_up, units = "days")) %>%
  filter(completion == 1) %>% 
  distinct()


library(lubridate)


base <- merge(base,
              base_temp,
              by = c("study","id","age","age_cat","sex","race","region"),
              all.x = T)#excludes all patients, which should not be analysed

base <- base %>%
  mutate(dis_duration = as.Date(start_date) - as.Date(dis_date)) %>%
  select(-dis_date)

# 1 = Caucasian
# 2 = Black
# 3 = Asian
# 7 = Native American
# 8 = Pacific islander
# 77 = Unknown
# 88 = Other

base$race <- as.factor(base$race)
levels(base$race) <- c("Caucasian","Black","Asian","Native American","Pacific islander","Other")
base$race <- as.character(base$race)
base$sex <- ifelse(base$sex == 1,"M","F")


#_____MRI

file_mri <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2309/Files/Analysis Ready Datasets/SAS_analysis/a_mri.sas7bdat"
CFTY720D2309_mri <- haven::read_sas(data_file = file_mri)

mri_temp <- data.frame(study = CFTY720D2309_mri$STY1A,
                       id = CFTY720D2309_mri$STYSID1A,
                       arm = CFTY720D2309_mri$TGPDSC1A,
                       para = CFTY720D2309_mri$MRIRSL1A,
                       unit = CFTY720D2309_mri$RSLUNT1A,
                       value = CFTY720D2309_mri$RSLVAL1N,
                       invalid_scan = CFTY720D2309_mri$FSTERPHA,
                       visit = CFTY720D2309_mri$VIS1N,
                       day = CFTY720D2309_mri$DAYNI_1N) %>%
  filter(invalid_scan == 0) %>%
  select(-invalid_scan)


# Take T2 lesion volume, T2 new or enlarged lesion count
# pivotize data set to wide format
mri_temp <- mri_temp %>%
  filter(para %in% c("T2_Vol","T2New_R")) %>%
  pivot_wider(names_from = c(para,unit),values_from = value) %>%
  select(-T2_Vol_)


names(mri_temp)[c(6,7)] <- c("t2_volume","t2_new_or_enlarged")
mri_temp <- mri_temp %>%
  mutate(unit = ifelse(!is.na(t2_volume),"mm?",NA))


mri_temp$id <- as.factor(mri_temp$id)

mri_temp <- merge(mri_temp,
                  base[,-which(names(base) == "n_pre_rel")],
                  by = c("study","id","arm"),
                  all.y = T) #excludes all patients, which should not be analysed


mri <- dplyr::bind_rows(mri,mri_temp) %>%
  as.data.frame()

#_____ARR  

file_arr <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2309/Files/Analysis Ready Datasets/SAS_analysis/a_rel.sas7bdat"
CFTY720D2309_arr <- haven::read_sas(data_file = file_arr)

arr_temp <- data.frame(study = CFTY720D2309_arr$STY1A,
                       id = CFTY720D2309_arr$STYSID1A,
                       arm = CFTY720D2309_arr$TGPDSC1A,
                       date_rel_start = CFTY720D2309_arr$RLPSTT1O,
                       date_rel_end = CFTY720D2309_arr$RLPEND1O,
                       relapse = CFTY720D2309_arr$RLPCFM1C,
                       comulated_rel = CFTY720D2309_arr$RELNUM2, #Number of confirmed relapses
                       day = CFTY720D2309_arr$DAY_1N,
                       follow_up = CFTY720D2309_arr$INSTUDY1, #days in study (non imputed)
                       arr_from_study = CFTY720D2309_arr$ARR1) # ARR of non imputed confirmed relapses

arr_temp <- merge(arr_temp,
                  base,
                  by = c("study","id","arm","follow_up"),
                  all.y = T)#excludes all patients, which should not be analysed

# Filter all relapses (relapse == 1) within clinicla trial (para == CLINICAL EPISODE)
# compute follow-up in days
# compute ARR
arr_temp <- arr_temp %>%
  #filter(analysis == "Y") %>% #exclude non Analysis/ITT population
  mutate(date_rel_start = as.Date(date_rel_start),
         date_rel_end = as.Date(date_rel_end)) %>%
  group_by(id) %>%
  mutate(comulated_rel = sum(relapse, na.rm = T),
         arr = comulated_rel/(follow_up/365.25)) %>%
  distinct() %>%
  filter(relapse > 0 | is.na(relapse)) # get rid of all patients without any relapse





arr <- bind_rows(arr,arr_temp)

#_____EDSS

file_edss <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2309/Files/Analysis Ready Datasets/SAS_analysis/a_edss.sas7bdat"
CFTY720D2309_edss <- haven::read_sas(data_file = file_edss)

edss_temp <- data.frame(study = CFTY720D2309_edss$STY1A,
                        id = CFTY720D2309_edss$STYSID1A,
                        arm = CFTY720D2309_edss$TGPDSC1A,
                        edss_score = CFTY720D2309_edss$EXDSC1N,
                        visual = CFTY720D2309_edss$FUNVSU1N,
                        brainstm = CFTY720D2309_edss$FUNBRA1N,
                        pyramidl = CFTY720D2309_edss$FUNPYR1N,
                        cerebell = CFTY720D2309_edss$FUNCER1N,
                        sensory = CFTY720D2309_edss$FUNSEN1N,
                        blad_bow = CFTY720D2309_edss$FUNBOW1N,
                        cerebral = CFTY720D2309_edss$FUNMEN1N,
                        walk = CFTY720D2309_edss$WLKDSN2C,
                        visit = CFTY720D2309_edss$VIS1N,
                        day = CFTY720D2309_edss$DAYNI_1N)


edss_temp$id <- as.factor(edss_temp$id)

edss_temp <- merge(edss_temp,
                  base[,-which(names(base) == "n_pre_rel")],
                  by = c("study","id","arm"),
                  all.y = T)#excludes all patients, which should not be analysed


edss <- dplyr::bind_rows(edss,edss_temp) %>%
  as.data.frame()




#_____Disease progression

file_dp <- "E:/SDT-RP-11223/files_NVT_SA_FTY720D2309/Files/Analysis Ready Datasets/SAS_analysis/a_msfc.sas7bdat"
CFTY720D2309_dp <- haven::read_sas(data_file = file_dp)

dp_temp <- data.frame(study = CFTY720D2309_dp$STY1A,
                      id = CFTY720D2309_dp$STYSID1A,
                      arm = CFTY720D2309_dp$TGPDSC1A,
                      peg_test = CFTY720D2309_dp$HPT,
                      x25_walk = CFTY720D2309_dp$T25FW,
                      visit = CFTY720D2309_dp$VIS1N,
                      day = CFTY720D2309_dp$DAYNI_1N) %>%
  unique()


dp_temp <- merge(dp_temp,
                 edss_temp %>% 
                   select(study,id,arm,edss_score,pyramidl,day,visit),
                 by = c("study","id","arm","day","visit"),
                 all.y = T)#excludes all patients, which should not be analysed


dp <- bind_rows(dp,dp_temp)



rm(dp_temp)
rm(CFTY720D2309_dp)
rm(file_dp)
rm(CFTY720D2309_mri)
rm(mri_temp)
rm(file_mri)
rm(CFTY720D2309_arr)
rm(arr_temp)
rm(file_arr)
rm(CFTY720D2309_edss)
rm(edss_temp)
rm(file_edss)
rm(file_base)
rm(CFTY720D2309_base)
rm(base)
rm(base_temp)
#####

#_______9.) WA21092
#####

#_____CREATION OF BASLIBE DATASET

file_base <- "E:/SDT-RP-11223/files_RCH-WA21092-PRIM1/Files/Analysis Ready Datasets/SAS_analysis/adsl.sas7bdat"
WA21092_base <- haven::read_sas(data_file = file_base)

base <- data.frame(   study = WA21092_base$STUDYID,
                      id = WA21092_base$USUBJID,
                      arm = WA21092_base$TRT01P,
                      age = WA21092_base$AGE,
                      sex = WA21092_base$SEX,
                      race = WA21092_base$RACE,
                      follow_up = WA21092_base$XTRTEDY,
                      analysis = WA21092_base$ITTOLEFL, #ITT population. Excludes small follow-ups and 
                      completion = WA21092_base$TRTCMPFL,
                      region = NA,
                      dis_duration = WA21092_base$ONSETDUR,
                      n_pre_rel = WA21092_base$RLP2YR) %>%
  distinct() %>%
  mutate(dis_duration = as.difftime(dis_duration, units = "days"),
         follow_up = as.difftime(follow_up, units = "days")) %>%
  filter(analysis == "Y")


#_____MRI

file_mri <- "E:/SDT-RP-11223/files_RCH-WA21092-PRIM2/Files/Analysis Ready Datasets/SAS_analysis/amri.sas7bdat"
WA21092_mri <- haven::read_sas(data_file = file_mri)

mri_temp <- data.frame(study = WA21092_mri$STUDYID,
                       id = WA21092_mri$USUBJID,
                       arm = WA21092_mri$TRT01P,
                       para = WA21092_mri$PARAM,
                       unit = WA21092_mri$AVALU,
                       value = WA21092_mri$AVAL,
                       visit = WA21092_mri$VISITNUM,
                       day = WA21092_mri$VISITDY)

unique(mri_temp$para)

# Take T2 lesion volume, T2 new or enlarged lesion count
# pivotize data set to wide format
 mri_temp <- mri_temp %>%
  filter(para %in% c( "T2 Lesion Count",
                      "T2 Lesion Volume",
                      "New/Enlarging T2 Lesion Count")) %>%
  pivot_wider(names_from = c(para,unit),values_from = value)

names(mri_temp)[6:8] <- c("t2_new_or_enlarged","t2_new_or_enlarged_bl","t2_volume")
mri_temp <- mri_temp %>%
  mutate(unit = ifelse(!is.na(t2_volume),"cm?",NA))

# Add baseline T2 Count to t2_new_or_enlarged:
mri_temp[which(mri_temp$visit == 2),"t2_new_or_enlarged"] <- mri_temp[which(mri_temp$visit == 2),"t2_new_or_enlarged_bl"]
mri_temp$t2_new_or_enlarged_bl <- NULL

mri_temp <- merge(mri_temp,
                  base %>% select(-n_pre_rel),
                  by = c("study","id","arm"),
                  all.y = T)#excludes all patients, which should not be analysed


mri <- mri %>%
  mutate(completion=ifelse(completion == 1,"Y","N"))

mri <- dplyr::bind_rows(mri,mri_temp) %>%
  as.data.frame()

#_____ARR  

file_arr <- "E:/SDT-RP-11223/files_RCH-WA21092-PRIM1/Files/Analysis Ready Datasets/SAS_analysis/aae.sas7bdat"
WA21092_arr <- haven::read_sas(data_file = file_arr)

arr_temp <- data.frame(study = WA21092_arr$STUDYID,
                       id = WA21092_arr$USUBJID,
                       arm = WA21092_arr$TRT01P,
                       follow_up = WA21092_arr$XTRTEDY,
                       #analysis_end_day = WA21092_arr$AENDY,
                       study_day_rel_start = WA21092_arr$XAESTDY,
                       study_day_rel_end = WA21092_arr$XAEEDY,
                       period = WA21092_arr$APERIODC,
                       relapse = WA21092_arr$AECAT,
                       day = WA21092_arr$XAESTDY) %>%
  distinct()

arr_temp <- merge(arr_temp %>% filter(relapse == "MS RELAPSE" & 
                                      period == "TREATMENT"), #Exclude observations after treatment phase (Exclusion of OLE (open lable extention)/SFU (safety) Phase)
                  base,
                  by = c("study","id","arm","follow_up"),
                  all.y = T) %>%#excludes all patients, which should not be analysed
  distinct()


# Filter all relapses (relapse == 1) within clinicla trial (para == CLINICAL EPISODE)
# compute follow-up in days
# compute ARR
arr_temp <- arr_temp %>%
  #filter(analysis == "Y") %>% # Exclude all non ITT/Analysis population
  #filter(!is.na(study_day_rel_start)) %>% # exclude observations with missing start days of relpase.
  #filter(!is.na(study_day_rel_end)) %>% # exclude observations with missing end days of relpase.
  group_by(id) %>%
  mutate(relapse = ifelse(relapse == "MS RELAPSE",1,relapse),
         comulated_rel = sum(relapse, na.rm = T),
         arr = comulated_rel/(follow_up/365.25))


arr$completion <- ifelse(arr$completion == 1,"Y","N")
arr <- bind_rows(arr,arr_temp)


#____Alternative way to get relapses:

#file_arr <- "E:/SDT-RP-11223/files_RCH-WA21092-PRIM2/Files/Analysis Ready Datasets/SAS_analysis/arlpit.sas7bdat"
#WA21092_arr <- haven::read_sas(data_file = file_arr)

#arr_temp <- data.frame(study = WA21092_arr$STUDYID,
#                       id = WA21092_arr$USUBJID,
#                       arm = WA21092_arr$TRT01P,
#                       study_day_rel_start = WA21092_arr$ADY,
#                       comulated_rel = WA21092_arr$NRLP, #relapse during treatment
#                       day = WA21092_arr$EXP96DYS) # Study day of end of relapse


# compute follow-up in days
# compute ARR
#arr_temp <- arr_temp %>%
#  group_by(id) %>%
#  mutate(follow_up = max(day),
#         arr = comulated_rel/(follow_up/365.25)) %>%
#  filter(comulated_rel > 0) # get rid of all patients without any relapse

#arr <- bind_rows(arr,arr_temp)


#_____EDSS

file_edss <- "E:/SDT-RP-11223/files_RCH-WA21092-PRIM1/Files/Analysis Ready Datasets/SAS_analysis/aedss.sas7bdat"
WA21092_edss <- haven::read_sas(data_file = file_edss)

edss_temp <- data.frame(study = WA21092_edss$STUDYID,
                        id = WA21092_edss$USUBJID,
                        arm = WA21092_edss$TRT01P,
                        para = WA21092_edss$PARAMCD,
                        value = WA21092_edss$AVAL,
                        visit = WA21092_edss$VISITNUM,
                        day = WA21092_edss$VISITDY)

unique(edss_temp$para)

edss_temp <- edss_temp %>%
  filter(!is.na(value)) %>%
  filter(!(is.na(day) & visit %in% c(44,999))) %>% # filter all unsceduld or "safety follow up" visits with no time point assigned
  filter(!(is.na(day) & is.na(visit))) %>% # filter all rows, with no day and visit entry (no assignment to timepoint possible)
  pivot_wider(names_from = para,values_from = value)


names(edss_temp)[6:14] <- c("walk",
                            "blad_bow",
                            "brainstm",
                            "cerebell",
                            "cerebral",
                            "edss_score",
                            "pyramidl",
                            "sensory",
                            "visual")

edss_temp <- merge(edss_temp,
                  base,
                  by = c("study","id","arm"),
                  all.y = T)#excludes all patients, which should not be analysed

edss$completion <- ifelse(edss$completion == 1,"Y","N")


edss <- bind_rows(edss,edss_temp)




#_____Disease progression

file_dp <- "E:/SDT-RP-11223/files_RCH-WA21092-PRIM1/Files/Analysis Ready Datasets/SAS_analysis/afcs.sas7bdat"
WA21092_dp <- haven::read_sas(data_file = file_dp)

dp_temp <- data.frame(study = WA21092_dp$STUDYID,
                      id = WA21092_dp$USUBJID,
                      arm = WA21092_dp$TRT01P,
                      para = WA21092_dp$PARAM,
                      value = WA21092_dp$AVAL,
                      visit = WA21092_dp$VISITNUM,
                      day = WA21092_dp$XADY,
                      period = WA21092_dp$APERIODC)

dp_temp <- dp_temp %>%
  filter(para %in% c("Average of 25 Foot Walk Tests","Average of 9 Hole Peg Tests")) %>%
  filter(period == "TREATMENT") %>%
  pivot_wider(names_from = para,values_from = value) %>%
  select(-period)

names(dp_temp)[c(6,7)] <- c("x25_walk","peg_test")

dp_temp <- merge(dp_temp,
                 edss_temp %>% 
                   select(study,id,arm,edss_score,pyramidl,day,visit),
                 by = c("study","id","arm","day","visit"),
                 all.y = T)#excludes all patients, which should not be analysed


dp <- bind_rows(dp,dp_temp)

rm(WA21092_dp)
rm(file_dp)
rm(dp_temp)
rm(WA21092_mri)
rm(mri_temp)
rm(file_mri)
rm(WA21092_arr)
rm(arr_temp)
rm(file_arr)
rm(WA21092_edss)
rm(edss_temp)
rm(file_edss)
rm(base)
rm(WA21092_base)
rm(file_base)
#####

#_______10.) WA21093
#####

#_____CREATION OF BASLIBE DATASET

file_base <- "E:/SDT-RP-11223/files_RCH-WA21093-PRIM1/Files/Analysis Ready Datasets/SAS_analysis/adsl.sas7bdat"
WA21093_base <- haven::read_sas(data_file = file_base)

base <- data.frame(   study = WA21093_base$STUDYID,
                      id = WA21093_base$USUBJID,
                      arm = WA21093_base$TRT01P,
                      age = WA21093_base$AGE,
                      sex = WA21093_base$SEX,
                      race = WA21093_base$RACE,
                      follow_up = WA21093_base$XTRTEDY,
                      analysis = WA21093_base$ITTOLEFL, #ITT population
                      completion = WA21093_base$TRTCMPFL,
                      region = NA,
                      dis_duration = WA21093_base$ONSETDUR,
                      n_pre_rel = WA21093_base$RLP2YR) %>%
  distinct() %>%
  mutate(dis_duration = as.difftime(dis_duration, units = "days"),
         follow_up = as.difftime(follow_up, units = "days")) %>%
  filter(completion == "Y")


#_____MRI

file_mri <- "E:/SDT-RP-11223/files_RCH-WA21093-PRIM2/Files/Analysis Ready Datasets/SAS_analysis/amri.sas7bdat"
WA21093_mri <- haven::read_sas(data_file = file_mri)

mri_temp <- data.frame(study = WA21093_mri$STUDYID,
                       id = WA21093_mri$USUBJID,
                       arm = WA21093_mri$TRT01P,
                       para = WA21093_mri$PARAM,
                       unit = WA21093_mri$AVALU,
                       value = WA21093_mri$AVAL,
                       visit = WA21093_mri$VISITNUM,
                       day = WA21093_mri$VISITDY)


# Take T2 lesion volume, T2 new or enlarged lesion count
# pivotize data set to wide format
mri_temp <- mri_temp %>%
  filter(!is.na(value)) %>%
  filter(visit != 999) %>% #visit 999 is an unsceduled visit. So, no time point can be assigned
  filter(para %in% c( "T2 Lesion Count",
                      "T2 Lesion Volume",
                      "New/Enlarging T2 Lesion Count")) %>%
  pivot_wider(names_from = c(para,unit),values_from = value) 

names(mri_temp)[6:8] <- c("t2_new_or_enlarged","t2_new_or_enlarged_bl","t2_volume")
mri_temp <- mri_temp %>%
  mutate(unit = ifelse(!is.na(t2_volume),"cm?",NA))

# Add baseline T2 Count to t2_new_or_enlarged:
mri_temp[which(mri_temp$visit == 2),"t2_new_or_enlarged"] <- mri_temp[which(mri_temp$visit == 2),"t2_new_or_enlarged_bl"]
mri_temp$t2_new_or_enlarged_bl <- NULL

mri_temp <- merge(mri_temp,
                   base,
                   by = c("study","id","arm"),
                   all.y = T)#excludes all patients, which should not be analysed


mri <- dplyr::bind_rows(mri,mri_temp) %>%
  as.data.frame()

#_____ARR  

file_arr <- "E:/SDT-RP-11223/files_RCH-WA21093-PRIM1/Files/Analysis Ready Datasets/SAS_analysis/aae.sas7bdat"
WA21093_arr <- haven::read_sas(data_file = file_arr)

arr_temp <- data.frame(study = WA21093_arr$STUDYID,
                       id = WA21093_arr$USUBJID,
                       arm = WA21093_arr$TRT01P,
                       follow_up = WA21093_arr$XTRTEDY,
                       #analysis_end_day = WA21093_arr$AENDY,
                       study_day_rel_start = WA21093_arr$XAESDY,
                       study_day_rel_end = WA21093_arr$XAEEDY,
                       period = WA21093_arr$APERIODC,
                       relapse = WA21093_arr$AECAT,
                       day = WA21093_arr$XAESDY) %>%
  distinct()

arr_temp <- merge(arr_temp %>% filter(relapse == "MS RELAPSE" & 
                                        period == "TREATMENT"), #Exclude observations after treatment phase (Exclusion of OLE (open lable extention)/SFU (safety) Phase)
                  base,
                  by = c("study","id","arm","follow_up"),
                  all.y = T) %>%#excludes all patients, which should not be analysed
  distinct()


# Filter all relapses (relapse == 1) within clinicla trial (para == CLINICAL EPISODE)
# compute follow-up in days
# compute ARR
arr_temp <- arr_temp %>%
  #filter(analysis == "Y") %>% # Exclude all non ITT/Analysis population
  #filter(!is.na(study_day_rel_start)) %>% # exclude observations with missing start days of relpase.
  #filter(!is.na(study_day_rel_end)) %>% # exclude observations with missing end days of relpase.
  group_by(id) %>%
  mutate(relapse = ifelse(relapse == "MS RELAPSE",1,relapse),
         comulated_rel = sum(relapse, na.rm = T),
         arr = comulated_rel/(follow_up/365.25))

arr <- bind_rows(arr,arr_temp)

#_____EDSS

file_edss <- "E:/SDT-RP-11223/files_RCH-WA21093-PRIM1/Files/Analysis Ready Datasets/SAS_analysis/aedss.sas7bdat"
WA21093_edss <- haven::read_sas(data_file = file_edss)

edss_temp <- data.frame(study = WA21093_edss$STUDYID,
                        id = WA21093_edss$USUBJID,
                        arm = WA21093_edss$TRT01P,
                        para = WA21093_edss$PARAMCD,
                        value = WA21093_edss$AVAL,
                        visit = WA21093_edss$VISITNUM,
                        day = WA21093_edss$VISITDY)


edss_temp <- edss_temp %>%
  filter(!is.na(value)) %>%
  filter(!(is.na(day) & visit %in% c(44,999))) %>% # filter all unsceduld or "safety follow up" visits with no time point assigned
  filter(!(is.na(day) & is.na(visit))) %>% # filter all rows, with no day and visit entry (no assignment to timepoint possible)
  pivot_wider(names_from = para,values_from = value)


names(edss_temp)[6:14] <- c("walk",
                            "blad_bow",
                            "brainstm",
                            "cerebell",
                            "cerebral",
                            "edss_score",
                            "pyramidl",
                            "sensory",
                            "visual")

edss_temp <- merge(edss_temp,
                  base,
                  by = c("study","id","arm"),
                  all.y = T)#excludes all patients, which should not be analysed

edss <- bind_rows(edss,edss_temp)



#_____Disease progression

file_dp <- "E:/SDT-RP-11223/files_RCH-WA21093-PRIM1/Files/Analysis Ready Datasets/SAS_analysis/afcs.sas7bdat"
WA21093_dp <- haven::read_sas(data_file = file_dp)

dp_temp <- data.frame(study = WA21093_dp$STUDYID,
                      id = WA21093_dp$USUBJID,
                      arm = WA21093_dp$TRT01P,
                      para = WA21093_dp$PARAM,
                      value = WA21093_dp$AVAL,
                      visit = WA21093_dp$VISITNUM,
                      day = WA21093_dp$XADY,
                      period = WA21093_dp$APERIODC)

dp_temp <- dp_temp %>%
  filter(para %in% c("Average of 25 Foot Walk Tests","Average of 9 Hole Peg Tests")) %>%
  filter(period == "TREATMENT") %>%
  pivot_wider(names_from = para,values_from = value) %>%
  select(-period)

names(dp_temp)[c(6,7)] <- c("x25_walk","peg_test")

dp_temp <- merge(dp_temp,
                 edss_temp %>% 
                   select(study,id,arm,edss_score,pyramidl,day,visit),
                 by = c("study","id","arm","day","visit"),
                 all.y = T)#excludes all patients, which should not be analysed


dp <- bind_rows(dp,dp_temp)

rm(WA21093_dp)
rm(file_dp)
rm(dp_temp)
rm(WA21093_mri)
rm(mri_temp)
rm(file_mri)
rm(WA21093_arr)
rm(arr_temp)
rm(file_arr)
rm(WA21093_edss)
rm(edss_temp)
rm(file_edss)
rm(file_base)
rm(WA21093_base)
rm(base)
#####

#_______11.) WA21493
#####

#_____CREATION OF BASLIBE DATASET

file_base <- "E:/SDT-RP-11223/files_RCH-WA21493-UPDT/Files/Analysis Ready Datasets/SAS_analysis/adsl.sas7bdat"
WA21493_base <- haven::read_sas(data_file = file_base)

base <- data.frame(   study = WA21493_base$STUDYID,
                      id = WA21493_base$USUBJID,
                      arm = WA21493_base$TRT01P,
                      age = WA21493_base$AGE,
                      sex = WA21493_base$SEX,
                      race = WA21493_base$RACE,
                      follow_up = WA21493_base$XTRTEDY,
                      analysis = WA21493_base$ITTOLEFL, #ITT population
                      completion = WA21493_base$TRTCMPFL,
                      region = NA,
                      dis_duration = WA21493_base$ONSETDUR,
                      n_pre_rel = WA21493_base$RLP2YR) %>%
  distinct() %>%
  mutate(dis_duration = as.difftime(dis_duration, units = "days"),
         follow_up = as.difftime(follow_up, units = "days")) %>%
  filter(completion == "Y")


#_____MRI

file_mri <- "E:/SDT-RP-11223/files_RCH-WA21493-UPDT/Files/Analysis Ready Datasets/SAS_analysis/amri.sas7bdat"
WA21493_mri <- haven::read_sas(data_file = file_mri)

mri_temp <- data.frame(study = WA21493_mri$STUDYID,
                       id = WA21493_mri$USUBJID,
                       arm = WA21493_mri$TRT01P,
                       para = WA21493_mri$PARAMCD,
                       unit = WA21493_mri$AVALU,
                       value = WA21493_mri$AVAL,
                       visit = WA21493_mri$VISITNUM,
                       period = WA21493_mri$APERIODC,
                       day = WA21493_mri$VISITDY) %>%
  filter(period == "TREATMENT") %>%
  select(-period)

unique(mri_temp$para)


# Take T2 lesion volume, T2 new or enlarged lesion count
# pivotize data set to wide format
mri_temp <- mri_temp %>%
  filter(visit != 999) %>% #visit 999 is an unsceduled visit. So, no time point can be assigned
  filter(para %in% c( "NEWT2",
                      "T2LESV",
                      "T2NEWENL",
                      "T2NEWLSN",
                      #"T2HYNUM", no observation occurs
                      "T2PERENL")) %>%
  filter(!is.na(value)) %>%
  pivot_wider(names_from = c(para,unit),values_from = c(value)) 

names(mri_temp)[6:10] <- c("t2_new_or_enlarged","t2_volume","t2_enlarged","t2_new","t2_per_enlarging")
mri_temp <- mri_temp %>%
  mutate(unit = ifelse(!is.na(t2_volume),"cm?",NA))




mri_temp <- merge(mri_temp,
                  base,
                  by = c("study","id","arm"),
                  all.y = T)#excludes all patients, which should not be analysed

mri <- dplyr::bind_rows(mri,mri_temp) %>%
  as.data.frame()

#_____ARR  

file_arr <- "E:/SDT-RP-11223/files_RCH-WA21493-UPDT/Files/Analysis Ready Datasets/SAS_analysis/aae.sas7bdat"
WA21493_arr <- haven::read_sas(data_file = file_arr)

arr_temp <- data.frame(study = WA21493_arr$STUDYID,
                       id = WA21493_arr$USUBJID,
                       arm = WA21493_arr$TRT01P,
                       follow_up = WA21493_arr$XTRTEDY,
                       follow_up_dis = WA21493_arr$TRTDISDY, #day of discontinuation in study
                       #analysis_end_day = WA21493_arr$AENDY,
                       study_day_rel_start = WA21493_arr$XAESTDY,
                       study_day_rel_end = WA21493_arr$XAEENDY,
                       period = WA21493_arr$APERIODC,
                       relapse = WA21493_arr$AECAT,
                       day = WA21493_arr$XAESTDY) %>%
  mutate(follow_up = ifelse(is.na(follow_up),follow_up_dis,follow_up)) %>% # replace NA in follow up by follow up of discontinuation.
  select(-follow_up_dis) %>%
  distinct()


# compute ARR
arr_temp <- merge(arr_temp %>% filter(relapse == "MS RELAPSE" & 
                                        period == "TREATMENT"), #Exclude observations after treatment phase (Exclusion of OLE (open lable extention)/SFU (safety) Phase)
                  base,
                  by = c("study","id","arm","follow_up"),
                  all.y = T) %>%#excludes all patients, which should not be analysed
  distinct()


# Filter all relapses (relapse == 1) within clinicla trial (para == CLINICAL EPISODE)
# compute follow-up in days
# compute ARR
arr_temp <- arr_temp %>%
  #filter(analysis == "Y") %>% # Exclude all non ITT/Analysis population
  #filter(!is.na(study_day_rel_start)) %>% # exclude observations with missing start days of relpase.
  #filter(!is.na(study_day_rel_end)) %>% # exclude observations with missing end days of relpase.
  group_by(id) %>%
  mutate(relapse = ifelse(relapse == "MS RELAPSE",1,relapse),
         comulated_rel = sum(relapse, na.rm = T),
         arr = comulated_rel/(follow_up/365.25))


arr <- bind_rows(arr,arr_temp)

#_____EDSS

file_edss <- "E:/SDT-RP-11223/files_RCH-WA21493-UPDT/Files/Analysis Ready Datasets/SAS_analysis/aedss.sas7bdat"
WA21493_edss <- haven::read_sas(data_file = file_edss)

edss_temp <- data.frame(study = WA21493_edss$STUDYID,
                        id = WA21493_edss$USUBJID,
                        arm = WA21493_edss$TRT01P,
                        para = WA21493_edss$PARAMCD,
                        value = WA21493_edss$AVAL,
                        visit = WA21493_edss$VISITNUM,
                        day = WA21493_edss$VISITDY)


edss_temp <- edss_temp %>%
  filter(!(is.na(day) & visit == 999)) %>% # filter all unsceduld visits with no time point assigned
  pivot_wider(names_from = para,values_from = value)


names(edss_temp)[6:13] <- c("blad_bow",
                            "brainstm",
                            "cerebell",
                            "cerebral",
                            "edss_score",
                            "pyramidl",
                            "sensory",
                            "visual")


edss_temp <- merge(edss_temp,
                  base,
                  by = c("study","id","arm"),
                  all.y = T)#excludes all patients, which should not be analysed


edss <- bind_rows(edss,edss_temp)

#_____Disease Progression

dp_temp <- edss_temp %>%
  select(study,id,arm,visit,day,edss_score,pyramidl)

dp <- bind_rows(dp,dp_temp)



rm(dp_temp)
rm(WA21493_mri)
rm(mri_temp)
rm(file_mri)
rm(WA21493_arr)
rm(arr_temp)
rm(file_arr)
rm(WA21493_edss)
rm(edss_temp)
rm(file_edss)
rm(base)
rm(file_base)
rm(WA21493_base)
#####

edss <- unique(edss)

# Adjust units of t2 volume

mri$unit <- ifelse(mri$unit == "cm?","cm3",mri$unit)
mri$unit <- ifelse(mri$unit == "mm?","mm3",mri$unit)





