setwd("C:/Users/Asus/Dropbox/CD_entropy/paper_con_template/esempi/MultivariateNormal")
results=read.csv("mvtnormal_MCMC.csv")
library(pdfCluster)

thinning=rep(c(1,0),25000)
thinning
which(thinning==1)
thinned_results=results[which(thinning==1),]
thinned_results
dd=pdfCluster::kepdf(thinned_results)

plot(density(thinned_results[,5]))
abline(v=0.95)
#maximum posterior in the null hypothesis set
idx=which((abs(dd@eval.points[,5]-0.9))<0.005)
max_=which.max(dd@estimate[idx])
max_
dd@eval.points[idx[max_],]
#value at maximum
dd@estimate[idx[max_]]
threshold=dd@estimate[idx[max_]]

new_post_=sample(1:nrow(dd@eval.points),prob = dd@estimate, replace = T, size = 100000)
new_post=cbind(dd@eval.points[new_post_,],dd@estimate[new_post_])
nrow(new_post) 

#sixth column contains the posterior value
1-mean(new_post[,6]>threshold)
 

#0.72




results_prof=read.csv("mvtnormal_profile_score.csv")
dd_prof=pdfCluster::kepdf(results_prof )
plot(dd_prof)
plot(density(results_prof[,1]))
dd_prof=density(results_prof[,1])
new_post_=sample(1:length(dd_prof$x),prob = dd_prof$y, replace = T, size = 10000)
new_post=cbind(dd_prof$x[new_post_],dd_prof$y[new_post_])
 
#sixth column contains the posterior value

idx=which((abs(dd_prof$x-0.9))<0.005)
max_=which.max(dd_prof$y[idx])
max_
dd_prof$x[idx[max_]]
#value at maximum
dd_prof$x[idx[max_]]
threshold=dd_prof$y[idx[max_]]
1-mean(new_post[,2]>threshold)
#0.08
plot(density(new_post[,1]))



library(ggplot2)
data <- data.frame(rho = c(results[,5], results_prof[,1]),
                   method=c(rep("marginal posterior", length(results[,5])),
                            rep("profiled matching", length(results_prof[,1]))))
# Create the Q-Q plot using ggplot2
p1=ggplot(data, aes(x=rho, color=method, legend=T)) +
  labs(x = " ", y = " ") +
  xlab(expression(rho))+
  geom_density(adjust=2)
p1  
ggtitle("Likelihood")



