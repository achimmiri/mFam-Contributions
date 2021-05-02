library(circlize)
library(devtools)
library(gdata)
library(jsonlite)
library(metaMS)
library(plotrix)
library(qdapRegex)
library(rJava)
library(rcdk)
library(rinchi)
library(rromeo)
library(readxl)
library(CAMERA)
library(RMassBank)
library(stringr)
library(stringi)
library(squash)
library(tools)
library(webchem)
library(zeallot)
library(BBmisc)
library(purrr)
library(schoolmath)
library(classyfireR)
library(sen2r)
########## Step1: Reading the API key
Sys.setenv(CHEMSPIDER_KEY = CHEMSPIDERKEY)
rr_auth(CHEMSPIDERKEY)
apikey=CHEMSPIDERKEY
#########################################################
#########################################################
args <- commandArgs(TRUE)
File1<-args[1]
##########################################################
Fi <- unlist(strsplit(File1, "/"))
Fi1 <-c(paste(Fi[-length(Fi)], collapse = "/"), last(Fi))
Fi2<-Fi1[1]
Fi3<-unlist(strsplit(Fi2, "/"))
Fi4<-c(paste(Fi3[-length(Fi3)], collapse = "/"), last(Fi3))
Fi5<-Fi4[1]
Fi6<-paste(Fi5,"raw data","exported as raw msp",sep="/")
### Input################
Fi7<-paste(Fi6,"/",sep="")
Fi8<-paste(Fi5,"converted to msp",sep="/")
### Output###############
Fi9<-paste(Fi8,"/",sep="")
#######################################################
####### Step2: Reading the All the required files
AIN<-read.table("ADI.txt",sep="\t",header=F,quote="",stringsAsFactors = FALSE)
####################Step3: READ the meta data file
RXF<-read_excel(File1, sheet = 1, col_names = TRUE,skip=1,.name_repair="minimal")
RXF[] <- lapply(RXF, function(x) type.convert(as.character(x)))
RXF3<- RXF
#######################################################
#Lmeda<-add(unname(table(as.character(unique(RXF3[["File"]])))))
Lmeda<-add(unname(table(as.character(RXF3[["File"]]))))
LmeCmu<-abs(cumsum(Lmeda))
LmeCmu1<-c(0,LmeCmu)
SFileNam<-as.character(unique(RXF3[["File"]]))
#######################################################
MaKComFile<-function(gDire)
{
  Fi <- unlist(strsplit(gDire, "/"))
  Fi1 <-c(paste(Fi[-length(Fi)], collapse = "/"), last(Fi))
  Fi2<-Fi1[1]
  Fi3<-unlist(strsplit(Fi2, "/"))
  Fi4<-c(paste(Fi3[-length(Fi3)], collapse = "/"), last(Fi3))
  Fi5<-Fi4[1]
  Fi6<-paste(Fi5,"raw data","exported as raw msp",sep="/")
  ### Input################
  Fi7<-paste(Fi6,"/",sep="")
  Fi8<-paste(Fi5,"converted to msp",sep="/")
  ### Output MSP file ... Input for combined
  Fi9<-paste(Fi8,"/",sep="")
  ##########################
  BNA<-basename(comsub(c(Fi6, Fi8), sep = ""))
  ### Combined outpuz
  Fout<-paste(Fi9,BNA,sep="")
  Fout1<-paste(Fout,"combined","msp",sep=".")
  
  ##########################
  if (!(length(list.files(path = Fi8, pattern = "\\.combined.msp$")) > 0))
  {
    LF<-Sys.glob(file.path(Fi9, "*.passed.msp"))
    myfiles = lapply(LF, readLines)
    lapply(myfiles, write, file=Fout1, append=T)
  }
}

########################################################
MaKlist<-function(gFile)
{
  
  lines <- readLines(gFile)
  lst <-split(lines, cumsum(lines==""))
  lst1 <-lapply(lst, function(x) if (x[1] == "") x[-1] else x)
  LL<-sapply(lst, length)
  IR<-which(unname(LL) == 0)
  if(length(IR)>0){LL1<-LL[-(IR)]}else{LL1 <- LL}
  ###### MSP List#############################
  lst2 <-lst1[names(LL1)]
  lst2<-purrr::compact(lst2)
  #################
  PMZL<-unname(sapply(lst2, function(x) grep("PRECURSORMZ",x)))
  PMZL1<-Filter(length,PMZL)
  ##### making addut information
  pattern <- "PRECURSORTYPE|ADDUCTIONNAME"
  PTYL<-unname(sapply(lst2, function(x) grep(pattern,x)))
  PTYL1<-Filter(length,PTYL)
  ##############################
  ###### Adduct list ###########
  faddu<-c()
  if(length(PTYL1)>0){
    for(i in c(1:length(PTYL1))){
      addu<-grep("PRECURSORTYPE:|ADDUCTIONNAME:",lst2[[i]]) 
      addu1<-stringr::str_trim(gsub('PRECURSORTYPE:|ADDUCTIONNAME:','',lst2[[i]][addu]))
      addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},error=function(cond){message("error happened in precursortype extraction")})
      addu3<-tryCatch({AIN[AIN$V1==addu2,]$V2},error=function(cond){message("adduct match is empty")})
      faddu<-c(faddu,addu3)
      
    }
  }else(print("AIn is empty"))
  #### The list with added adduct value 
  lst3<-tryCatch({Map(c, lst2, faddu)},error=function(cond){message("faddu is empty")})
  ######### Precursor MZ information
  lst4<-c()
  for(i in 1:length(lst3)){
    LV<-tryCatch({lst3[[i]]},error=function(cond){message("list is empty")})
    LaV<-tryCatch({BBmisc::getLast(LV)},error=function(cond){message("last element is empty")})
    SR<-tryCatch({stringi::stri_startswith_fixed(LV, 'PRECURSORMZ:')},error=function(cond){message("PRECURSORMZ is empty")})
    ##SR<-tryCatch({stri_detect_fixed(LV,"PRECURSORMZ:")},error=function(cond){message("PRECURSORMZ is empty")})
    SR1<-tryCatch({which(SR)},error=function(cond){message("index is empty")})
    SR2<-tryCatch({LV[[SR1]]},error=function(cond){message("index does not exist")})
    PV<-tryCatch({stringr::str_trim(stringr::str_replace(SR2, "PRECURSORMZ:", ""))},error=function(cond){message("index does not exist")})
    #S1r<-tryCatch({str_replace(LaV, "M", PV)},error=function(cond){message("replacement does not exist")})
    S1r<-tryCatch({stringr::str_replace(LaV, "M", PV)},error=function(cond){message("replacement does not exist")})
    S2r<-tryCatch({as.numeric(pander::evals(S1r)[[1]]$result)},error=function(cond){message("result is empty")})
    lst4<-c(lst4,S2r)
  }
  ######################################################
  #### adding the adduct adjusted mass to make final list
  lst5<-tryCatch({Map(c, lst3, lst4)},error=function(cond){message("Either lst3 or lst4 is empty")})
  ######## Precursormz value 
  fmass<-c()
  for(i in c(1:length(PMZL1))){
    LV1<-tryCatch({lst3[[i]]},error=function(cond){message("list is empty")})
    #S1R<-tryCatch({stri_detect_fixed(LV1,"PRECURSORMZ:")},error=function(cond){message("PRECURSORMZ is empty")})
    #stringi::stri_startswith_fixed(LV, 'PRECURSORMZ:')
    S1R<-tryCatch({stringi::stri_startswith_fixed(LV1, 'PRECURSORMZ:') },error=function(cond){message("PRECURSORMZ is empty")})
    S1R1<-tryCatch({which(S1R)},error=function(cond){message("Index is empty")})
    S1R2<-tryCatch({LV1[S1R1]},error=function(cond){message("Element is empty")})
    PV1<-tryCatch({stringr::str_trim(stringr::str_replace(S1R2, "PRECURSORMZ:", ""))},error=function(cond){message("PRECURSORMZ empty")})
    PV2<-tryCatch({as.numeric(PV1)},error=function(cond){message("PRECURSORMZ does not exists")})
    fmass<-c(fmass,PV2)
  }
  ####################################################
  lst7<-tryCatch({Map(c, lst5, fmass)},error=function(cond){message("Either lst5 or fmass is empty")})
  ######### Coping to new list
  lst8<-lst7 
  ##### Adduct adjusted mass values
  AAMZV<-sapply(lst5, function(x) BBmisc::getLast(x))
  AAMZV <- setNames(as.numeric(AAMZV), names(AAMZV))
  ###############################
  ###### Retention time list#####
  ### Making the RT list######### 
  RTL<-sapply(lst2, function(x) grep("RETENTIONTIME",x))
  IND<-unname(RTL)
  ############# Final Retention Test list
  FRTL <- c()
  if(length(IND)>0){
    for(i in c(1:length(IND))){
      FRTL <- c(FRTL,get(names(lst2)[i],lst2)[IND[i]])
    }}
  ##############################################
  FRTL1<-tryCatch({as.numeric(sapply(strsplit(FRTL, ":"),`[`, 2))},error = function(cond){message("out is empty")})
  ##############################################
  return(list(lst2, PMZL1, PTYL1,faddu,lst3,lst4,lst5,fmass,lst7,lst8,AAMZV,RTL,IND,FRTL,FRTL1))
  
}
#########################################################
#NFFilter(RRV,AIN,lst2,fmass,FRTL1)
NFFilter<-function(InMEDA,InAdVA,InMSPL,InPMZ,InRTL)  
{
  ##############
  out<-c()
  ##############
  IV<-as.character(InMEDA[["InChI"]])
  EM<-as.character(InMEDA["Exact mass"])
  EM1<-as.numeric(EM)
  ADV<-as.character(InMEDA[["Adduct"]])
  ADV1<-tryCatch({qdapRegex::ex_between(ADV, "[", "]")[[1]]},error=function(cond){message("adduct value is missing")})
  ADV2<-tryCatch({InAdVA[InAdVA$V1==ADV1,]$V2},error=function(cond){message("adduct value matching is not found")})
  ADV3<-tryCatch({stringr::str_replace(ADV2, "M",as.character(EM1))},error=function(cond){message("adduct value replacement is not found")})
  ADV4<-tryCatch({as.numeric(pander::evals(ADV3)[[1]]$result)},error=function(cond){message("getting the result")})
  #####################
  PPm=ADV4*(10/(1000000))
  ####################
  MPPmL=ADV4-PPm
  MPPmU=ADV4+PPm
  ####################
  
  Tmass<-InPMZ[InPMZ >= MPPmL & InPMZ <= MPPmU]
  ITmass<-which(InPMZ %in% Tmass)
  #####################
  VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
  VRTL<-VRT-0.08
  VRTU<-VRT+0.08
  ######################
  
  TRTL<-InRTL[InRTL >=VRTL & InRTL <= VRTU]
  ITRTL<-which(InRTL %in% TRTL)
  ######################
  INLL<-intersect(ITmass,ITRTL)
  ########################################
  if(length(ITRTL) >= 1){
    if(length(ITmass) >= 1){
      if(length(INLL) == 1){
        #####################
        
        F1FPL<-InMSPL[INLL]
        F2FPL<-F1FPL
        #######################
        Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("there is an error in list")})
        ########################
        FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
        ########################
        FNA1<-which(stri_detect_fixed(FNA,"NAME:"))
        FNA2<-InMEDA[["Name"]]
        FNA3<-as.character(FNA2)
        FNAM<-paste("NAME:",FNA3,sep=" ")
        out<-c(out,FNAM)
        ####################
        FRA1<-which(stri_detect_fixed(FNA,"RETENTIONTIME:"))
        F1RA1<-FNA[FRA1]
        out<-c(out,F1RA1)
        ##################
        FMZ1<-which(stri_detect_fixed(FNA,"PRECURSORMZ:"))
        F1MZ1<-FNA[FMZ1]
        out<-c(out,F1MZ1)
        ###################
        FPT1<-which(stri_detect_fixed(FNA,"PRECURSORTYPE:"))
        F1PT1<-FNA[FPT1]
        out<-c(out,F1PT1)
        ###########################
        FIN1<-InMEDA[["Ionization mode"]]
        F1IN1<-as.character(FIN1)
        F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
        out<-c(out,F2IN1)
        ###########################
        F1ONT<-paste("Ontology:",as.character(InMEDA[["Compound class"]]),sep=" ")
        out<-c(out,F1ONT)
        ###########################
        if(!is.na(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI="))
        {
          FINK<-paste("INCHIKEY:","",sep=" ")
          out<-c(out,FINK)
          FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
          out<-c(out,FINCH)
        }else if(!is.na(as.character(InMEDA[["InChI"]]))){
          
          if(tryCatch({webchem::is.inchikey(as.character(InMEDA[["InChI"]]))},error=function(cond){message("webchem is not able to fecth information")}))
          {
            FINK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
            out<-c(out,FINK)
            FINCH<-paste("INCHI:","",sep=" ")
            out<-c(out,FINCH)
          }
        }else{
          FINK<-paste("INCHIKEY:","",sep=" ")
          FINCH<-paste("INCHI:","",sep=" ")
          out<-c(out,FINK)
          out<-c(out,FINCH)
        }
        ######################################
        FSIM<-paste("SMILES:",as.character(InMEDA[["SMILES"]]),sep=" ")
        out<-c(out,FSIM)
        ##############################
        FFOR<-InMEDA[["Formula"]]
        FFOR1<-paste("FORMULA:",FFOR,sep=" ")
        out<-c(out,FFOR1)
        ############################
        FINS<-which(stri_detect_fixed(FNA,"INTENSITY:"))
        FINS1<-FNA[FINS]
        out<-c(out,FINS1)
        ###########################
        FAUT<-as.character(InMEDA[["Authors"]])
        FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
        out<-c(out,FAUT1)
        #############################
        FLIC<-paste("LICENSE:",sep=" ")
        out<-c(out,FLIC)
        ###########################
        FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
        out<-c(out,FCIE)
        ##########################
        FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
        FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
        out<-c(out,FINST1)
        #########################
        FINS<-as.character(InMEDA[["INSTRUMENT"]])
        FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
        out<-c(out,FINS1)
        #########################
        FCOM<-paste("COMMENT:")
        out<-c(out,FCOM)
        ########################
        FNPA<-which(stri_detect_fixed(FNA,"Num Peaks:"))
        F1NPA<-FNA[FNPA]
        out<-c(out,F1NPA)
        #######################
        Fpea<-FNA[(FNPA+1):Find]
        ######################
        if(is.na(Fpea))
        {
          Fpea1<-FNA[(FNPA+1)]
          out<-c(out,Fpea1)
          
        }else{
          MV=ADV4
          tes1<-unlist(strsplit(Fpea, "\t"))
          tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
          tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
          tes4<-which(tes2 > (3+MV))
          if(length(tes4)>1)
          {
            tes5<-tes2[-tes4]
            tes6<-tes3[-tes4]
            tes7<-paste(tes5,tes6,sep="\t")
            out<-c(out,tes7)
          }else{
            out<-c(out,Fpea)
          }
        }
        
      ##########################################
      }else if(length(INLL) > 1){
        ##########
        MONMS=InMSPL[INLL]
        TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
        TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
        TRA2<-abs(VRT-TRA1)
        TRA3<-which.min(TRA2)
        TRA4<-INLL[TRA3]
        TRA5<-InMSPL[TRA4]
        ###########################
        F1FPL<-TRA5
        F2FPL<-F1FPL
        ##########################
        Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("there is an error in list")})
        ########################
        FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
        ########################
        FNA1<-which(stri_detect_fixed(FNA,"NAME:"))
        FNA2<-InMEDA[["Name"]]
        FNA3<-as.character(FNA2)
        FNAM<-paste("NAME:",FNA3,sep=" ")
        out<-c(out,FNAM)
        ####################
        FRA1<-which(stri_detect_fixed(FNA,"RETENTIONTIME:"))
        F1RA1<-FNA[FRA1]
        out<-c(out,F1RA1)
        ##################
        FMZ1<-which(stri_detect_fixed(FNA,"PRECURSORMZ:"))
        F1MZ1<-FNA[FMZ1]
        out<-c(out,F1MZ1)
        ###################
        FPT1<-which(stri_detect_fixed(FNA,"PRECURSORTYPE:"))
        F1PT1<-FNA[FPT1]
        out<-c(out,F1PT1)
        ###########################
        FIN1<-InMEDA[["Ionization mode"]]
        F1IN1<-as.character(FIN1)
        F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
        out<-c(out,F2IN1)
        ###########################
        F1ONT<-paste("Ontology:",as.character(InMEDA[["Compound class"]]),sep=" ")
        out<-c(out,F1ONT)
        ###########################
        if(!is.na(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI="))
        {
          FINK<-paste("INCHIKEY:","",sep=" ")
          out<-c(out,FINK)
          FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
          out<-c(out,FINCH)
        }else if(!is.na(as.character(InMEDA[["InChI"]]))){
          
          if(tryCatch({webchem::is.inchikey(as.character(InMEDA[["InChI"]]))},error=function(cond){message("webchem is not able to fecth information")}))
          {
            FINK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
            out<-c(out,FINK)
            FINCH<-paste("INCHI:","",sep=" ")
            out<-c(out,FINCH)
          }
        }else{
          FINK<-paste("INCHIKEY:","",sep=" ")
          FINCH<-paste("INCHI:","",sep=" ")
          out<-c(out,FINK)
          out<-c(out,FINCH)
        }
        ###############################################
        FSIM<-paste("SMILES:",as.character(InMEDA[["SMILES"]]),sep=" ")
        out<-c(out,FSIM)
        ############################################
        FFOR<-InMEDA[["Formula"]]
        FFOR1<-paste("FORMULA:",FFOR,sep=" ")
        out<-c(out,FFOR1)
        ############################################
        FINS<-which(stri_detect_fixed(FNA,"INTENSITY:"))
        FINS1<-FNA[FINS]
        out<-c(out,FINS1)
        ###########################################
        FAUT<-as.character(InMEDA[["Authors"]])
        FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
        out<-c(out,FAUT1)
        ###########################################
        FLIC<-paste("LICENSE:",sep=" ")
        out<-c(out,FLIC)
        ##########################################
        FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
        out<-c(out,FCIE)
        #########################################
        FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
        FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
        out<-c(out,FINST1)
        #########################################
        FINS<-as.character(InMEDA[["INSTRUMENT"]])
        FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
        out<-c(out,FINS1)
        ######################
        FCOM<-paste("COMMENT:")
        out<-c(out,FCOM)
        ####################
        FNPA<-which(stri_detect_fixed(FNA,"Num Peaks:"))
        F1NPA<-FNA[FNPA]
        out<-c(out,F1NPA)
        ###################
        Fpea<-FNA[(FNPA+1):Find]
        ###################
        if(is.na(Fpea))
        {
          Fpea1<-FNA[(FNPA+1)]
          out<-c(out,Fpea1)
          
        }else{
          MV=ADV4
          tes1<-unlist(strsplit(Fpea, "\t"))
          tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
          tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
          tes4<-which(tes2 > (3+MV))
          if(length(tes4)>1)
          {
            tes5<-tes2[-tes4]
            tes6<-tes3[-tes4]
            tes7<-paste(tes5,tes6,sep="\t")
            out<-c(out,tes7)
          }else{
            out<-c(out,Fpea)
          }
        } ### this is the else MV= ADV4...closing 
      ###################  
      } else {
        #print("entering the line 458")
        PASS<-RRV
      }  ## main else part
      ###testing if this works
      return(out)
      #####################
    } #ITmass
  } #ITRTL
} 

##################### Function InchiKey filter
##############################################
#Ikfilter(IV,lst2,RRV,AIN,fmass,FRTL1)
Ikfilter <- function(InKeyVal,InMSPL,InMEDA,InAdVA,InPMZ,InRTL){
  #################
  out<-c()
  #################
  IINF<-tryCatch({webchem::cts_compinfo(InKeyVal)},error=function(cond){message("webchecm could not fetch the info")})
  ##################
  if(length(IINF)> 0){
    ###################
    IK<-tryCatch({IINF[[1]][1]},error=function(cond){message("Inchiley value is empty")})
    PMZ<-tryCatch({IINF[[1]][4]},error=function(cond){message("PrecursorMZ value is empty")})
    FM<-tryCatch({IINF[[1]][5]},error=function(cond){message("Formula value is empty")})
    CID<-tryCatch({webchem::get_cid(IK, from = "inchikey")},error=function(cond){message("webchecm could not fetch the info")})
    CID1<-tryCatch({CID%>% dplyr::select(cid)},error=function(cond){message("CompoundID is empty; check previous step")})
    CID2<-as.character(CID1)
    CID3<-gsub("[[:punct:]]", "",CID2 )
    CID4<-unlist(strsplit(CID3, " "))
    CID5<-paste(CID4, collapse = ';')
    ###########################
    SM<-tryCatch({webchem::cir_query(IK,"smiles")},error=function(cond){message("webchecm could not fetch the info")})
    SM1<-tryCatch({SM[[1]][1]},error=function(cond){message("smiles Information fetch error")})
    ###########################
    InchiV<-tryCatch({rinchi::get.inchi(SM1)},error=function(cond){message("webchecm could not fetch the info")})
    ############################
    AUIN<-as.character(InMEDA[["Adduct"]])
    ############################
    if(length(AUIN) > 0 & !is.na(PMZ) ){
      ##########################
      AUIN1<-tryCatch({qdapRegex::ex_between(AUIN, "[", "]")[[1]]},error=function(cond){message("Adduct value is missing")})
      AUIN2<-tryCatch({InAdVA[InAdVA$V1==AUIN1,]$V2},error=function(cond){message("Adduct value is missing")})
      AAMS<-tryCatch({stringr::str_replace(AUIN2, "M",as.character(PMZ$exactmass))},error=function(cond){message("Missing adduct replacement")})
      AAMS1<-tryCatch({as.numeric(pander::evals(AAMS)[[1]]$result)},error=function(cond){message("Error in adduct replacement step")})
      #########################
      PPm=AAMS1*(10/(1000000))
      #########################
      MPPmL=AAMS1-PPm
      MPPmU=AAMS1+PPm
      ###########################
      Tmass<-InPMZ[InPMZ >= MPPmL & InPMZ <= MPPmU]
      ITmass<-which(InPMZ %in% Tmass)
      ###########################
      VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
      VRTL<-VRT-0.08 
      VRTU<-VRT+0.08
      #############################
      TRTL<-InRTL[InRTL >=VRTL & InRTL <= VRTU]
      ITRTL<-which(InRTL %in% TRTL)
      ################################
      if(length(ITRTL) >= 1){
        if(length(ITmass) >= 1){
          ##################
          INLL<-intersect(ITmass,ITRTL)
          ###################
          if(length(INLL) == 1){
            ###################
            F1FPL<-InMSPL[INLL]
            #####################
            SM1<-as.character(InMEDA[["SMILES"]])
            #####################
            F2FPL<-tryCatch({F1FPL},error=function(cond){message("List value is empty")})
            ######################
            Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("List value is empty")})
            ########################
            FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
            ########################
            FNA1<-which(stri_detect_fixed(FNA,"NAME:"))
            FNA2<-InMEDA[["Name"]]
            FNA3<-as.character(FNA2)
            #######################
            FNAM<-paste("NAME:",FNA3,sep=" ")
            out<-c(out,FNAM)
            ########################
            FRA1<-which(stri_detect_fixed(FNA,"RETENTIONTIME:"))
            F1RA1<-FNA[FRA1]
            out<-c(out,F1RA1)
            ################################
            FMZ1<-which(stri_detect_fixed(FNA,"PRECURSORMZ:"))
            F1MZ1<-FNA[FMZ1]
            out<-c(out,F1MZ1)
            ##########################
            FPT1<-which(stri_detect_fixed(FNA,"PRECURSORTYPE:"))
            F1PT1<-FNA[FPT1]
            out<-c(out,F1PT1)
            ###########################
            FIN1<-InMEDA[["Ionization mode"]]
            F1IN1<-as.character(FIN1)
            F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
            out<-c(out,F2IN1)
            #########################
            IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},error=function(cond){message("Classifier could not fecth the information")})
            ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},error=function(cond){message("Classifier could not fecth the information")})), sep = ","))
            F1ONT<-paste("Ontology:",ONTV,sep=" ")
            out<-c(out,F1ONT)
            ##############################
            FINK<-paste("INCHIKEY:",tryCatch({IK$inchikey},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
            out<-c(out,FINK)
            FINCH<-paste("INCHI:",InchiV,sep=" ")
            out<-c(out,FINCH)
            FSIM<-paste("SMILES:",SM1,sep=" ")
            out<-c(out,FSIM)
            ##################
            FFOR<-FM$formula
            FFOR1<-paste("FORMULA:",FFOR,sep=" ")
            out<-c(out,FFOR1)
            #######################
            FINS<-which(stri_detect_fixed(FNA,"INTENSITY:"))
            FINS1<-FNA[FINS]
            out<-c(out,FINS1)
            #######################
            FAUT<-as.character(InMEDA[["Authors"]])
            FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
            out<-c(out,FAUT1)
            ##################
            FLIC<-paste("LICENSE:",sep=" ")
            out<-c(out,FLIC)
            ##################
            FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
            out<-c(out,FCIE)
            ##################
            FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
            FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
            out<-c(out,FINST1)
            #####################
            FINS<-as.character(InMEDA[["INSTRUMENT"]])            
            FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
            out<-c(out,FINS1)
            ####################
            FCOM<-paste("COMMENT:")
            out<-c(out,FCOM)
            ##################
            FNPA<-which(stri_detect_fixed(FNA,"Num Peaks:"))
            F1NPA<-FNA[FNPA]
            out<-c(out,F1NPA)
            ###################
            Fpea<-FNA[(FNPA+1):Find]
            #########################
            if(is.na(Fpea))
            {
              Fpea1<-FNA[(FNPA+1)]
              
              
            }else{
              
              MV=AAMS1
              tes1<-unlist(strsplit(Fpea, "\t"))
              tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
              tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
              tes4<-which(tes2 > (3+MV))
              if(length(tes4)>1)
              {
                tes5<-tes2[-tes4]
                tes6<-tes3[-tes4]
                tes7<-paste(tes5,tes6,sep="\t")
                out<-c(out,tes7)
              }else{
                out<-c(out,Fpea)
                
              }
            }
            
          ####################################  
          } else if(length(INLL) > 1){                          
            #print("entering the 496")
            MONMS=InMSPL[INLL]
            TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
            TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
            TRA2<-abs(VRT-TRA1)
            TRA3<-which.min(TRA2)
            TRA4<-INLL[TRA3]
            TRA5<-InMSPL[TRA4]
            #######################
            F1FPL<-TRA5
            #####################
            SM1<-as.character(InMEDA[["SMILES"]])
            #InMEDA[["SMILES"]]<-SM1
            #InMEDA[["PubChem CID"]]<-CID5
            #####################
            F2FPL<-F1FPL
            ######################
            Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("List value is empty")})
            ########################
            FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("List value is empty")})
            ########################
            FNA1<-which(stri_detect_fixed(FNA,"NAME:"))
            FNA2<-InMEDA[["Name"]]
            FNA3<-as.character(FNA2)
            ######################
            FNAM<-paste("NAME:",FNA3,sep=" ")
            out<-c(out,FNAM)
            ######################
            FRA1<-which(stri_detect_fixed(FNA,"RETENTIONTIME:"))
            F1RA1<-FNA[FRA1]
            out<-c(out,F1RA1)
            ################################
            FMZ1<-which(stri_detect_fixed(FNA,"PRECURSORMZ:"))
            F1MZ1<-FNA[FMZ1]
            out<-c(out,F1MZ1)
            ################################
            FPT1<-which(stri_detect_fixed(FNA,"PRECURSORTYPE:"))
            F1PT1<-FNA[FPT1]
            out<-c(out,F1PT1)
            #############################
            FIN1<-InMEDA[["Ionization mode"]]
            F1IN1<-as.character(FIN1)
            F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
            out<-c(out,F2IN1)
            ###############################
            IKCRV<-classyfireR::get_classification(InKeyVal)
            ONTV<-do.call(paste, c(as.list(IKCRV@classification$Classification), sep = ","))
            F1ONT<-paste("Ontology:",ONTV,sep=" ")
            out<-c(out,F1ONT)
            ##############################
            FINK<-paste("INCHIKEY:",tryCatch({IK$inchikey},error=function(cond){message("Inchikey value is empty")}),sep=" ")
            out<-c(out,FINK)
            FINCH<-paste("INCHI:",InchiV,sep=" ")
            out<-c(out,FINCH)
            FSIM<-paste("SMILES:",SM1,sep=" ")
            out<-c(out,FSIM)
            ##############################
            FFOR<-FM$formula
            FFOR1<-paste("FORMULA:",FFOR,sep=" ")
            out<-c(out,FFOR1)
            ###############################
            FINS<-which(stri_detect_fixed(FNA,"INTENSITY:"))
            FINS1<-FNA[FINS]
            out<-c(out,FINS1)
            #############################
            FAUT<-as.character(InMEDA[["Authors"]])
            FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
            out<-c(out,FAUT1)
            #############################
            FLIC<-paste("LICENSE:",sep=" ")
            out<-c(out,FLIC)
            #############################
            FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
            out<-c(out,FCIE)
            ############################
            FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
            FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
            out<-c(out,FINST1)
            ##########################
            FINS<-as.character(InMEDA[["INSTRUMENT"]])
            FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
            out<-c(out,FINS1)
            ########################
            FCOM<-paste("COMMENT:")
            out<-c(out,FCOM)
            #######################
            FNPA<-which(stri_detect_fixed(FNA,"Num Peaks:"))
            F1NPA<-FNA[FNPA]
            out<-c(out,F1NPA)
            #######################
            Fpea<-FNA[(FNPA+1):Find]
            ########################
            if(is.na(Fpea))
            {
              Fpea1<-FNA[(FNPA+1)]
              
              
            }else{
              
              MV=AAMS1
              tes1<-unlist(strsplit(Fpea, "\t"))
              tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
              tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
              tes4<-which(tes2 > (3+MV))
              if(length(tes4)>1)
              {
                tes5<-tes2[-tes4]
                tes6<-tes3[-tes4]
                tes7<-paste(tes5,tes6,sep="\t")
                out<-c(out,tes7)
              }else{
                out<-c(out,Fpea)
                
              }
            }
          ################################  
          } else{
            #print("entering the line 750")
            PASS<-RRV
          }
          ####### This is testing , if this works
          return(out)
          #################
        } ### this is mz closing brace
      } ## this is RT closing braces
    }
  }
}
##############################################
##for(i in 1:1)
for(i in 1:length(LmeCmu1))
{
  
  Val=LmeCmu1[i]
  Val1=LmeCmu1[i+1]
  nVal=Val+1
  if(!is.na(Val) & !is.na(Val1))
  {
    #print("entering the line ...781")
    NRXF3<-RXF3[nVal:Val1,]
    FiNA<-SFileNam[i]
    BNFiNA<-basename(FiNA)
    BNFiNAEX<-tools::file_path_sans_ext(BNFiNA)
    FINAMSP<-paste(BNFiNAEX,"msp",sep=".")
    ###############FiNA1 ... INput ..msp file
    #####OUNA2 ... Output file
    FiNA1<-list.files(path =Fi6 , pattern = FINAMSP, recursive = TRUE, full.names = TRUE)
    ###############
    OUNA<-tools::file_path_sans_ext(FiNA)
    OUNA1<-paste(OUNA,"passed","msp",sep = ".")
    OUNA2<-paste(Fi9,OUNA1,sep="")
    #########################
    if((length(FiNA1) >= 1))
    {
      #print("entering the line ...line 797")
      #####################################
      CaIF<-tryCatch({MaKlist(FiNA1[1])},error=function(cond){message("some mistake happened in filename..file name must be empty")})
      lst2<-tryCatch({CaIF[[1]]},error=function(cond){message("some mistake happened in indexes")})
      fmass<-tryCatch({CaIF[[8]]},error=function(cond){message("some mistake happened in indexes")})
      FRTL1<-tryCatch({CaIF[[15]]},error=function(cond){message("some mistake happened in indexes")})
      ########################################
      if(file.exists(FiNA1[1]))
      {
       # print("entering the line ...806")
        ######################################
        if(length(NRXF3)>0){
          for (i in 1:length(NRXF3)){
            ##############
            RRV<-NRXF3[i,]
            ################
            PIN<-ifelse(is.na(RRV[["InChI"]]), 1, 0)
            PSM<-ifelse(is.na(RRV[["SMILES"]]), 1, 0)
            PPC<-ifelse(is.na(RRV[["PubChem CID"]]), 1, 0)
            ####################
            if(PIN == 1 & PSM == 1 & PPC == 1){
              #########
              PASS<-RRV
              #########
            }
            else{
              #print("enter the line ...839")
              if(!is.na(as.character(RRV[["InChI"]])) & !startsWith(as.character(RRV[["InChI"]]),'not available') & !startsWith(as.character(RRV[["InChI"]]),'CAS:')){
                #print("enter the line ...841")
                if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(RRV[["InChI"]])))},error=function(cond){message("some mistake happened in filename..file name must be empty")})){
                  IV<-stringr::str_trim(as.character(RRV[["InChI"]]))
                  if(!is.na(IV)){
                  outn<-Ikfilter(IV,lst2,RRV,AIN,fmass,FRTL1)
                  len1<-length(outn)
                  ############################
                  if(len1 > 1)
                  {
                    #print("enter the line ...845")
                    cat(sapply(outn, toString), file=OUNA2, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                    
                  }else{
                    FINKE<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                    ### adding this code new 
                    len3<-length(FINKE)
                    if(len3 > 1){
                      #print("enter the line ...856")
                      cat(sapply(FINKE, toString), file=OUNA2, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                    } ###len3 is end
                    
                  } # end of the else loop
                  
                  } # end of is.na(IV)
                } # end of is.inchikey
                  else{
                  #print("enter the line ...868")
                  FINKE1<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                  ### adding this code new 
                  len3<-length(FINKE1)
                  if(len3 > 1){
                    #print("enter the line ...873")
                    cat(sapply(FINKE1, toString), file=OUNA2, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                  } 
                }
                
              ####################################
              }else if(!is.na(as.character(RRV[["InChI"]])) & startsWith(as.character(RRV[["InChI"]]),'InChI=')){
                #print("enter the line ...861")
                IV1<-stringr::str_trim(as.character(RRV[["InChI"]]))
                FSMV<-tryCatch({chemspiderapi::post_convert(IV1,inputFormat = "InChI",outputFormat ="SMILES", apikey <- apikey)},error=function(cond){message("webchecm could not fetch the info")})
                FSMV1<-unname(FSMV)
                FSMV2<-tryCatch({rinchi::get.inchi.key(FSMV1)},error=function(cond){message("webchecm could not fetch the info")})
                ##########
                if(!is.na(FSMV2)){
                outn1<-Ikfilter(FSMV2,lst2,RRV,AIN,fmass,FRTL1)
                len2<-length(outn1)
                #########################
                if(len2 > 1)
                {
                  #print("enter the line ...878")
                  cat(sapply(outn1, toString), file=OUNA2, sep="\n",append=TRUE)
                  cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                  
                } # end of if loop
                } # FSMV2 is empty
                else{
                  FINKE1<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                  ### adding this code new 
                  len3<-length(FINKE1)
                  if(len3 > 1){
                    #print("enter the line ...888")
                    cat(sapply(FINKE1, toString), file=OUNA2, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                  } ###len3 is end
                } # end of the else loop
              ##############################
              }else if(!is.na(as.character(RRV[["InChI"]])) & startsWith(as.character(RRV[["InChI"]]),'CAS:')){
                #print("enter the line ...884")
                CV<-as.character(RRV[["InChI"]])
                CV1<-stringr::str_replace(CV,pattern='CAS:',replacement ="")
                CV2<-stringr::str_trim(as.character(CV1))
                ####################
                FIV=tryCatch({aw_query(CV2, from = 'cas')[[1]]},error=function(cond){message("webchecm could not fetch the info")})
                ##############################
                if(!is.na(FIV))
                {
                  FIV1<-FIV$inchikey
                  if(!is.na(FIV1)){
                    outn3<-Ikfilter(FIV1,lst2,RRV,AIN,fmass,FRTL1)
                    len4<-length(outn3)
                    if(len4 > 1)
                    {
                      #print("enter the line ...912")
                      cat(sapply(outn3, toString), file=OUNA2, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                      
                    } # end of len4
                  } ## FIV is empty
                    else{
                      FSMIL<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                      len3<-length(FSMIL)
                      if(len3 > 1){
                       # print("enter the line ...921")
                        cat(sapply(FSMIL, toString), file=OUNA2, sep="\n",append=TRUE)
                        cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                      }
                    } # end of else loop
                    #########
                #################################
                }else{
                  
                  FCAS<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                  ### adding this code new 
                  len4<-length(FCAS)
                  if(len4 > 1){
                    #print("enter the line ...935")
                    cat(sapply(FCAS, toString), file=OUNA2, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                  }
                }
              #################################
              }else if(!is.na(as.character(RRV[["SMILES"]])) & !startsWith(as.character(RRV[["SMILES"]]),'not available')){
                #print("enter the line ...915")
                F1SM<-as.character(RRV[["SMILES"]])
                F1SM1<-tryCatch({rinchi::get.inchi.key(F1SM)},error=function(cond){message("webchecm could not fetch the info")})
                #############
                ### adding this to F1SM1 is not empty
                if(!is.na(F1SM1)){
                  outn4<-Ikfilter(F1SM1,lst2,RRV,AIN,fmass,FRTL1)
                  loutn5<-length(outn4)
                  ################
                  if(loutn5 > 1)
                  {
                    #print("enter the line ...954")
                    cat(sapply(outn4, toString), file=OUNA2, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                    
                  }else{
                    #print("enter the line ...959")
                    FSMIL<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                    len3<-length(FSMIL)
                    if(len3 > 1){
                      #print("enter the line ...963")
                      cat(sapply(FSMIL, toString), file=OUNA2, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                    } ### end of len3
                  } ## end of else loop
                  ##### added this new else loop
                } else{
                  #print("enter the line ...970")
                  FSMIL<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                  len3<-length(FSMIL)
                  if(len3 > 1){
                    #print("enter the line ...974")
                    cat(sapply(FSMIL, toString), file=OUNA2, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                  } ### end of len3
                }
                
              ###################################
              }else if(!is.na(as.character(RRV[["PubChem CID"]])) & !startsWith(as.character(RRV[["PubChem CID"]]),'not available')){
                #print("enter the line ...939")
                FPUCID<-as.character(RRV[["PubChem CID"]])
                FINSM<-tryCatch({webchem::cs_convert(as.numeric(FPUCID),from = "csid", to = "smiles")},error=function(cond){message("webchecm could not fetch the info")})
                ###########################
                if(!is.na(FINSM))
                {
                  FIINK<-tryCatch({rinchi::get.inchi.key(FINSM)},error=function(cond){message("webchecm could not fetch the info")})
                  if(!is.na(FIINK)){
                    outn6<-Ikfilter(FIINK,lst2,RRV,AIN,fmass,FRTL1)
                    len7<-length(outn6)
                    if(len7 > 1)
                    {
                      #print("enter the line ...994")
                      cat(sapply(outn6, toString), file=OUNA2, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                      
                    }else{
                      #print("enter the line ...999")
                      FCIDR<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                      len3<-length(FCIDR)
                      if(len3 > 1){
                        #print("enter the line ...1003")
                        cat(sapply(FCIDR, toString), file=OUNA2, sep="\n",append=TRUE)
                        cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                      } ## end of len3
                    } # end of else loop
                    
                  }#FIINK is empty..inchikey is empty
                  else{
                    #print("enter the line ...1011")
                    FSMIL<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                    len3<-length(FSMIL)
                    if(len3 > 1){
                      #print("enter the line ...1014")
                      cat(sapply(FSMIL, toString), file=OUNA2, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                    } ### end of len3
                  }# end of else loop
                  #############################
                }else{
                  #print("enter the line ...1037")
                  FSMIL<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                  len3<-length(FSMIL)
                  if(len3 > 1){
                    #print("enter the line ...1014")
                    cat(sapply(FSMIL, toString), file=OUNA2, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA2, sep="\n",append=TRUE)
                  } ### end of len3
                } 
              }
              ###############################
              else{
                PASS1 <-RRV
              }
                
              }
            } ### big else loop
          } ## for loop NRXF3
        } ## End ...NRXF3
        
      } ## FiNA is closing one
      ###################################
    }
  }
