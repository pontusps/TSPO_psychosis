

            ############################################ 
# ##########          Systematic Pubmed search          ########## #
            ############################################

            # Pontus P. Sigray, KI, Stockholm, February 2017

#Load packages
library(RISmed)
library(xlsx)

#Specify search 
search_topic.pet <- "((( (psychosis) OR (psychotic disorder) OR (schizophrenia) )) AND ((translocator protein) OR (translocator protein 18 kDa) OR (TSPO) OR (TSPO18kDA) OR (peripheral benzodiazepine receptor) OR (PBR))) AND ((positron emission tomography) OR (PET))"

#Specify search dates
startdate<-'2004/01/01' #year when the first 2nd generation TSPO tracer was initially published on: http://onlinelibrary.wiley.com/doi/10.1002/syn.20027/abstract
stopdate<-'2017/02/20'

#Perform search
search_query <- EUtilsSummary(search_topic.pet, mindate=startdate, maxdate=stopdate)

#Retrieve pubmed's search syntax
search_terms<-QueryTranslation(search_query)

#Retrieve pmid of all hits
pmids<-QueryId(search_query)
n.hits<-length(pmids)

#Retrieve first and last author of all hits
records<- EUtilsGet(search_query)
author.list<-Author(records)

first.athor<-sapply(author.list,function(x) x$LastName[1])
last.author<-sapply(author.list,function(x) x$LastName[nrow(x)])

#Create output data sheet
pubmed_data<-data.frame('Title'=ArticleTitle(records),'Abstract'=AbstractText(records),First.Author=first.athor,Last.Author=last.author,pmids=PMID(records))

pubmed_data$Title<-as.character(pubmed_data$Title)
pubmed_data$Abstract<-as.character(pubmed_data$Abstract)

#Write retrieved data from pubmed to an excel file
write.xlsx(pubmed_data,file = paste0("./DerivedData/Systematic_literature_search/pubmed_dat.xlsx"),row.names = F)

#Write pubmed search terms to .txt file.  
write(as.character(search_terms),file = paste0("./DerivedData/Systematic_literature_search/pubmed_search_terms.txt"),ncolumns = 1)
