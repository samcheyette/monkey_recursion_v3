library(ggplot2)
library(reshape)
library(grid)
library(dplyr)

##########READ FILE#########################################
draw_from <- ""
inp_file <- (gsub(" ", "", paste(draw_from, "output2.csv")))

data = read.csv(inp_file)
m = melt(data, id=c("who", "hyp_dist", "HorM", "which", "P",
			 "alpha", "corr"))


t1 <- 	theme(axis.text=element_text(size=18), 
		strip.text.x = element_text(size = 20),
		plot.title=element_text(size=18),
		axis.text.x=element_blank(),
		#axis.text.x=element_text(size=13, angle=45),
		axis.title.x=element_blank(),
		axis.text.y=element_text(size=17),

		legend.title=element_blank(),
		legend.text=element_text(size=17),
		 legend.key.size = unit(4, 'lines'))


############################################################
###################PLOT1####################################

m.1 <- m %>% 
		filter(grepl("mod_trce", as.character(HorM))) %>%

		mutate(variable=(factor(gsub("OT", 
			"OTHER", as.character(variable))))) %>%
		mutate(variable=(factor(gsub("CE", 
			"CENTER-EMBEDDING", as.character(variable))))) %>%
		mutate(variable=(factor(gsub("CR", 
			"CROSSING", as.character(variable))))) %>%
		mutate(variable=(factor(gsub("TL", 
			"OTHER", as.character(variable)))))

		#mutate(variable=(factor(gsub("TL", 
			#"TAIL-RECURSION", as.character(variable)))))  




m.1 <- m.1 %>% 


			mutate(value=P*value) %>%
			ungroup %>%
			group_by(who) %>%
			#10 %>%
			top_n(20, wt=value) %>% 
			mutate(var2=paste(which, 
				as.character(as.numeric(P)))) %>%
			#group_by(Participant) %>%
			transform(var2 = factor(reorder(var2, -P))) %>%

			ungroup
		
head(m.1)

p.1 <- ggplot(data=m.1, aes(x=var2, y=value, group=variable)) +
		geom_bar(stat='identity', aes(fill=variable)) +
			geom_text(aes(label=which,x=var2, y=0.0),
				 size=8.0, angle=45) +
		facet_wrap(~who, scales="free_x")

p.1 <- p.1 + t1 + ylab("Probability")

ggsave("hyps.pdf", width=20, height=15) 


###################################################################
##########################DISTR####################################



m.dist <- m %>% 
			filter(grepl("dist", as.factor(hyp_dist))) %>%

			mutate(variance= P * (1 - P) ) %>%
			mutate(upper= (P + variance) * as.numeric(grepl("hum", HorM))) %>% 
			
			mutate(lower= (P - variance) * as.numeric(grepl("hum", HorM))) 





m.dist <- m.dist %>% 
		mutate(P2=(P * (grepl("hu",
			 as.character(HorM)))))


m.dist <- m.dist %>% 
		mutate(P3=(P * (grepl("mod",
			 as.character(HorM)))))


#m.dist <- m.dist %>% group_by(who, which) %>%
				#	mutate(psum=sum(P3)) 

					 #%>% 
					#filter(psum > thresh)
m.dist <- m.dist %>% group_by(who, which) %>%
					mutate(psum=sum(P2)) 

					 #%>% 
					#filter(psum > thresh)

m.dist <- m.dist %>% group_by(who, HorM) %>%

					top_n(30, wt=psum) %>%
					mutate(which2=paste(which, 
					as.character(as.numeric(psum)))) %>%
					transform(which2=factor(reorder(which2,-psum)))
head(m.dist)


p.dist.1 <- ggplot(m.dist, aes(x=which2, y=P, 
				group=HorM)) +
			geom_bar(stat='identity',position='dodge',
				aes(group=HorM, fill=HorM)) +
			#geom_errorbar(stat='identity', position='dodge',
					#aes(x=which2, ymin=lower, ymax=upper), 
					#size=0.5, alpha=0.9) +
			geom_text(aes(label=which,x=which2,y=0),
				 size=8.0, angle=45, vjust=-0.5) +
			facet_wrap(~who,
						scales="free_x")

p.dist.1 <- p.dist.1 +  
				theme(axis.text.x =
				 element_text(angle = 90, hjust = -1)) +
				t1  + ylab("Probability")


ggsave("dist1.pdf", width=20, height=12)
