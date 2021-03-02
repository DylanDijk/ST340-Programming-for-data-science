library(rio)
load('data/20newsgroups.rdata')

# documents matrix with (i,j)th entry indicator for the presence of the jth word in the ith post
# newsgroups vector    ith entry denotes tru value label for the ith post (which of the four topics it belongs to)
# groupnames names of the four topics
# wordlist list of the 100 words

groupnames
wordlist[100]



source('data/em_mixture_bernoullis.R')
source('data/my bernoulli EM function.R')

sz = 3000

smpl = sample(nrow(newsgroups), size=sz)
xs = documents[smpl,]
nrow(xs)

out = em_mix_bernoulli(xs, K = 4)

sum(abs(out$gammas - newsgroups.onehot[smpl,]))/sz



xs = documents
out_Q1bi = em_mix_bernoulli(xs, K  = 4)



for(i in 1:24)
sum(abs(out_Q1bi$gammas - newsgroups.onehot))/nrow(newsgroups)


library(gtools)
perm = permutations(n = 4, r = 4, v = 1:4)
gamma_minus_label = vector(length = 24)
for(i in 1:24){
  gamma_minus_label[i] = sum(abs(out_Q1bi$gammas[,perm[i,]] - newsgroups.onehot))
}

algorithm_labels = perm[which.min(gamma_minus_label),]

# This the average probability that the algorithm assigns to the "correct" label
sum((out_Q1bi$gammas[,algorithm_labels])[cbind(seq_len(nrow(newsgroups)), newsgroups)])/nrow(newsgroups)

library(doBy)
wordlist[which.maxn(out_Q1bi$mus[1,], n = 8)]
wordlist[which.maxn(out_Q1bi$mus[2,], n = 8)]
wordlist[which.maxn(out_Q1bi$mus[3,], n = 8)]
wordlist[which.maxn(out_Q1bi$mus[4,], n = 8)]









xss = matrix(c(1,0,0,1,1,0), nrow =3, ncol = 2, byrow=T)
em_mix_bernoulli(xss, K = 2)
my_em_mix_bernoulli(xss, K = 3)


out = my_em_mix_bernoulli(xs, K = 4)






head(out$gammas, 30)
head(newsgroups.onehot[smpl,], 30)
