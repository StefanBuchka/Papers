
#_____creation of common visit variable

#_____MRI
#####
mri$visit_com <- NA

visit_tmp <- mri[which(mri$study == "CFTY720D1201"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(0,3,4,5,6,777,999)
visit_tmp <- as.numeric(as.character(visit_tmp))
mri[which(mri$study == "CFTY720D1201"),"visit_com"] <- visit_tmp

visit_tmp <- mri[which(mri$study == "CFTY720D1201E"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(0,3,4,5,9,12,777,778,779,999)
visit_tmp <- as.numeric(as.character(visit_tmp))
mri[which(mri$study == "CFTY720D1201E"),"visit_com"] <- visit_tmp


visit_tmp <- mri[which(mri$study == "CFTY720D2301"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(0,6,12,777,777,778)
visit_tmp <- as.numeric(as.character(visit_tmp))
mri[which(mri$study == "CFTY720D2301"),"visit_com"] <- visit_tmp


visit_tmp <- mri[which(mri$study == "CFTY720D2302"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(0,777,777,777,999)
visit_tmp <- as.numeric(as.character(visit_tmp))
mri[which(mri$study == "CFTY720D2302"),"visit_com"] <- visit_tmp


visit_tmp <- mri[which(mri$study == "CFTY720D2309"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(0,6,12,501,6,12,24,777,800)
visit_tmp <- as.numeric(as.character(visit_tmp))
mri[which(mri$study == "CFTY720D2309"),"visit_com"] <- visit_tmp


visit_tmp <- mri[which(mri$study == "CAMMS223"),"visit"]
mri[which(mri$study == "CAMMS223"),"visit_com"] <- visit_tmp


visit_tmp <- mri[which(mri$study == "CAMMS323"),"visit"]
mri[which(mri$study == "CAMMS323"),"visit_com"] <- visit_tmp


visit_tmp <- mri[which(mri$study == "CAMMS324"),"visit"]
mri[which(mri$study == "CAMMS324"),"visit_com"] <- visit_tmp


visit_tmp <- mri[which(mri$study == "WA21092"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(0,6,11,22,33,433,999)
visit_tmp <- as.numeric(as.character(visit_tmp))
mri[which(mri$study == "WA21092"),"visit_com"] <- visit_tmp


visit_tmp <- mri[which(mri$study == "WA21093"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(0,6,11,22,33,433)
visit_tmp <- as.numeric(as.character(visit_tmp))
mri[which(mri$study == "WA21093"),"visit_com"] <- visit_tmp


visit_tmp <- mri[which(mri$study == "WA21493"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(0,1,2,3,4,5,6,22) #set basline to screening. Otherwise no baseline would be availlable
visit_tmp <- as.numeric(as.character(visit_tmp))
mri[which(mri$study == "WA21493"),"visit_com"] <- visit_tmp
#####



#_____EDSS and Diseas Progression
#####

edss$visit_com <- NA


visit_tmp <- edss[which(edss$study == "CAMMS223"),"visit"]
edss[which(edss$study == "CAMMS223"),"visit_com"] <- visit_tmp

visit_tmp <- edss[which(edss$study == "CAMMS323"),"visit"]
edss[which(edss$study == "CAMMS323"),"visit_com"] <- visit_tmp

visit_tmp <- edss[which(edss$study == "CAMMS324"),"visit"]
edss[which(edss$study == "CAMMS324"),"visit_com"] <- visit_tmp

visit_tmp <- edss[which(edss$study == "CFTY720D1201"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(0,3,6,777,999) #set screening to baseline measurement. Otherwise no baseline EDSS available
visit_tmp <- as.numeric(as.character(visit_tmp))
edss[which(edss$study == "CFTY720D1201"),"visit_com"] <- visit_tmp

visit_tmp <- edss[which(edss$study == "CFTY720D1201E"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(0,3,4,6,9,12,777,778,779,999) #set screening to baseline measurement. Otherwise no baseline EDSS available
visit_tmp <- as.numeric(as.character(visit_tmp))
edss[which(edss$study == "CFTY720D1201E"),"visit_com"] <- visit_tmp


visit_tmp <- edss[which(edss$study == "CFTY720D2301"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(1111,0,1,3,6,9,12,15,18,21,501,777,778,999)
visit_tmp <- as.numeric(as.character(visit_tmp))
edss[which(edss$study == "CFTY720D2301"),"visit_com"] <- visit_tmp


visit_tmp <- edss[which(edss$study == "CFTY720D2302"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(1111,0,3,6,9,501,601,602,603,777,999)
visit_tmp <- as.numeric(as.character(visit_tmp))
edss[which(edss$study == "CFTY720D2302"),"visit_com"] <- visit_tmp


visit_tmp <- edss[which(edss$study == "CFTY720D2309"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(1111,0,3,6,9,12,15,18,21,501,601,602,603,604,605,606,607,777,999)
visit_tmp <- as.numeric(as.character(visit_tmp))
edss[which(edss$study == "CFTY720D2309"),"visit_com"] <- visit_tmp


visit_tmp <- edss[which(edss$study == "WA21092"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(1111,0,3,6,9,11,14,17,19,22,100,103,109,111,117,433)
visit_tmp <- as.numeric(as.character(visit_tmp))
edss[which(edss$study == "WA21092"),"visit_com"] <- visit_tmp

visit_tmp <- edss[which(edss$study == "WA21093"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(1111,0,3,6,9,11,14,17,19,22,100,103,109,111,117,433)
visit_tmp <- as.numeric(as.character(visit_tmp))
edss[which(edss$study == "WA21093"),"visit_com"] <- visit_tmp

visit_tmp <- edss[which(edss$study == "WA21493"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(1111,0,3,6,9,11,14,17,19,22,25,28,
                       236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,
                       556,557,134,139,144,150,155,161,166,1066)
visit_tmp <- as.numeric(as.character(visit_tmp))
edss[which(edss$study == "WA21493"),"visit_com"] <- visit_tmp




dp$visit_com <- NA

visit_tmp <- dp[which(dp$study == "CAMMS223"),"visit"]
dp[which(dp$study == "CAMMS223"),"visit_com"] <- visit_tmp

visit_tmp <- dp[which(dp$study == "CAMMS323"),"visit"]
dp[which(dp$study == "CAMMS323"),"visit_com"] <- visit_tmp

visit_tmp <- dp[which(dp$study == "CAMMS324"),"visit"]
dp[which(dp$study == "CAMMS324"),"visit_com"] <- visit_tmp

visit_tmp <- dp[which(dp$study == "CFTY720D1201"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(1111,3,6,777,999)
visit_tmp <- as.numeric(as.character(visit_tmp))
dp[which(dp$study == "CFTY720D1201"),"visit_com"] <- visit_tmp

visit_tmp <- dp[which(dp$study == "CFTY720D1201E"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(1111,3,9,12,777,778,779,999)
visit_tmp <- as.numeric(as.character(visit_tmp))
dp[which(dp$study == "CFTY720D1201E"),"visit_com"] <- visit_tmp


visit_tmp <- dp[which(dp$study == "CFTY720D2301"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(1111,0,1,3,6,9,12,15,18,21,501,777,778,999)
visit_tmp <- as.numeric(as.character(visit_tmp))
dp[which(dp$study == "CFTY720D2301"),"visit_com"] <- visit_tmp


visit_tmp <- dp[which(dp$study == "CFTY720D2302"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(1111,0,3,6,9,501,601,602,603,777,999)
visit_tmp <- as.numeric(as.character(visit_tmp))
dp[which(dp$study == "CFTY720D2302"),"visit_com"] <- visit_tmp


visit_tmp <- dp[which(dp$study == "CFTY720D2309"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(1111,0,3,6,9,12,15,18,21,501,601,602,603,604,605,606,607,777,999)
visit_tmp <- as.numeric(as.character(visit_tmp))
dp[which(dp$study == "CFTY720D2309"),"visit_com"] <- visit_tmp


visit_tmp <- dp[which(dp$study == "WA21092"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(1111,0,3,6,9,11,14,17,19,22,100,103,106,109,117,433)
visit_tmp <- as.numeric(as.character(visit_tmp))
dp[which(dp$study == "WA21092"),"visit_com"] <- visit_tmp

visit_tmp <- dp[which(dp$study == "WA21093"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(1111,0,3,6,9,11,14,17,19,22,100,103,109,117,123,433)
visit_tmp <- as.numeric(as.character(visit_tmp))
dp[which(dp$study == "WA21093"),"visit_com"] <- visit_tmp

visit_tmp <- dp[which(dp$study == "WA21493"),"visit"]
visit_tmp <- as.factor(visit_tmp)
levels(visit_tmp) <- c(0,1,3,6,9,11,14,17,19,22,25,28,
                       236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,
                       556,557,134,139,144,150,155,161,166,1066)
visit_tmp <- as.numeric(as.character(visit_tmp))
dp[which(dp$study == "WA21493"),"visit_com"] <- visit_tmp
#####

rm(visit_tmp)
